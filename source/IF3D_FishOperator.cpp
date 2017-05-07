//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_FishOperator.h"
#include "IF3D_FishLibrary.h"
#include <chrono>
IF3D_FishOperator::IF3D_FishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_ObstacleOperator(grid, parser), theta_internal(0), angvel_internal(0),
sim_time(0), sim_dt(0), adjTh(0), adjDy(0), myFish(nullptr), angvel_integral{0,0,0},
new_curv(0), old_curv(0), new_Tp(0), ptrUinf_copy(nullptr)
{
	volume=0;
	for(int i=0;i<3;i++) transVel[i]=0;
	for(int i=0;i<3;i++) angVel[i]=0;
	for(int i=0;i<6;i++) J[i]=0;
}

IF3D_FishOperator::~IF3D_FishOperator()
{
	if(myFish not_eq nullptr) delete myFish;
}

void IF3D_FishOperator::create(const int step_id,const Real time, const Real dt, const Real *Uinf)
{
	// STRATEGY
	// we need some things already
	// - the internal angle at the previous timestep, obtained from integrating the actual def velocities
	// 						 (not the imposed deformation velocies, because they dont have zero ang mom)
	// - the internal angular velocity at previous timestep

	// 1. create midline
	// 2. integrate to find CoM, angular velocity, etc
	// 3. shift midline to CoM frame: zero internal linear momentum and angular momentum

	// 4. split the fish into segments (according to s)
	// 5. rotate the segments to computational frame (comp CoM and angle)
	// 6. for each Block in the domain, find those segments that intersect it
	// 7. for each of those blocks, allocate an ObstacleBlock

	// 8. put the 3D shape on the grid: SDF-P2M for sdf, normal P2M for udef
	// 9. create the Chi out of the SDF. In same sweep, compute the actual CoM
	// 10. compute all shit: linear momentum, angular momentum etc.
	// 11. correct deformation velocity to nullify momenta for the final discrete representation
	std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1, t23, t45, t67, t8, t910, t11;

	//MPI_Barrier(grid->getCartComm());
	t0 = std::chrono::high_resolution_clock::now();
	//const int Nsegments = NPPSEG;
	const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
	const Real dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
	const Real target_ds = vInfo[0].h_gridpoint/TGTPPB;
	const Real target_Nm = length/target_ds;
	////const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG) + 1;
	//const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1)) + 1;
	//const int Nsegments = (Nm-1)/(Nextension+1);
	/*
		- VolumeSegment_OBB's volume cannot be zero
		- therefore no VolumeSegment_OBB can be only occupied by extension midline
			points (which have width and height = 0)
		- performance of create seems to decrease if VolumeSegment_OBB are bigger
		- this is the smallest number of VolumeSegment_OBB (Nsegments) and points in
			the midline (Nm) to ensure at least one non ext. point inside all segments
	 */
	const int Nsegments = std::ceil(target_Nm/(Nextension+1));
	const int Nm = (Nextension+1)*Nsegments + 1;
	assert((Nm-1)%Nsegments==0);

	if (bCorrectTrajectory && time>0.312)
	{
		assert(followX < 0 && followY < 0);
		const Real velx_tot = Uinf[0] - transVel[0];
		const Real vely_tot = Uinf[1] - transVel[1];
		const Real AngDiff  = std::atan2(vely_tot,velx_tot);
		adjTh = (1.-dt) * adjTh + dt * AngDiff;
		const Real INST = (AngDiff*angVel[2]>0) ? AngDiff*std::fabs(angVel[2]) : 0;
		const Real PID = 0.1*adjTh + 0.01*INST;
		myFish->_correctTrajectory(PID, time, dt);
	}

	if (followX > 0 && followY > 0) //then i control the position
	{
		assert(not bCorrectTrajectory);
		const Real velx_tot = Uinf[0] - transVel[0];
		const Real vely_tot = Uinf[1] - transVel[1];
		const Real AngDiff  = std::atan2(vely_tot,velx_tot);

		// Control posDiffs
		const Real xDiff = (position[0] - followX)/length;
		const Real yDiff = (position[1] - followY)/length;
		
		adjTh = (1.-dt) * adjTh + dt * AngDiff;
		adjDy = (1.-dt) * adjDy + dt * yDiff;
		
		//If angle is positive: positive curvature only if Dy<0 (must go up)
		//If angle is negative: negative curvature only if Dy>0 (must go down)
		const Real INST = (AngDiff*yDiff<0) ? AngDiff*std::fabs(yDiff) : 0;
		const Real PID = 1*(adjTh-adjDy) + 1.*INST;
		myFish->_correctTrajectory(PID, time, dt);

		// Linearly increase (or decrease) amplitude to 1.2X (decrease to 0.8X)
		//(experiments observed 1.2X increase in amplitude when swimming faster)
		//if fish falls back 1 body length. Beyond that, will still increase but dunno if will work
		const Real amplitudeFactor = 1.*xDiff + 1.0;
		myFish->_correctAmplitude(amplitudeFactor, time, dt);
	}

	// 1.
	myFish->computeMidline(time);
	#ifdef __useSkin_
	myFish->computeSurface();
	#endif
	//MPI_Barrier(grid->getCartComm());
	t1 = std::chrono::high_resolution_clock::now();

	// 2. & 3.
	{
		volume_internal = myFish->integrateLinearMomentum(CoM_internal, vCoM_internal);
		assert(volume_internal > std::numeric_limits<Real>::epsilon());
		myFish->changeToCoMFrameLinear(CoM_internal, vCoM_internal);

		angvel_internal_prev = angvel_internal;
		myFish->integrateAngularMomentum(angvel_internal);
		J_internal = myFish->J;
		// update theta now with new angvel info
		//theta_internal -= 0.5*sim_dt*(angvel_internal+angvel_internal_prev);//negative: we subtracted this angvel
		myFish->changeToCoMFrameAngular(theta_internal, angvel_internal);
	}
	#ifndef NDEBUG
	{
		Real dummy_CoM_internal[2], dummy_vCoM_internal[2], dummy_angvel_internal;
		// check that things are zero
		const Real volume_internal_check = myFish->integrateLinearMomentum(dummy_CoM_internal,dummy_vCoM_internal);
		myFish->integrateAngularMomentum(dummy_angvel_internal);

		assert(std::abs(dummy_CoM_internal[0])<10*std::numeric_limits<Real>::epsilon());
		assert(std::abs(dummy_CoM_internal[1])<10*std::numeric_limits<Real>::epsilon());
		assert(std::abs(myFish->linMom[0])<10*std::numeric_limits<Real>::epsilon());
		assert(std::abs(myFish->linMom[1])<10*std::numeric_limits<Real>::epsilon());
		assert(std::abs(myFish->angMom)<10*std::numeric_limits<Real>::epsilon());
		assert(std::abs(volume_internal - volume_internal_check) < 10*std::numeric_limits<Real>::epsilon());
	}
	#endif
	//MPI_Barrier(grid->getCartComm());
	#ifdef __useSkin_
	myFish->surfaceToCOMFrame(theta_internal,CoM_internal);
	#endif
	t23 = std::chrono::high_resolution_clock::now();


	// 4. & 5.
	std::vector<VolumeSegment_OBB> vSegments(Nsegments);
	#pragma omp parallel for
	for(int i=0;i<Nsegments;++i) {
		const int next_idx = (i+1)*(Nm-1)/Nsegments;
		const int idx = i * (Nm-1)/Nsegments;
		// find bounding box based on this
		Real bbox[3][2] = {{1e9, -1e9}, {1e9, -1e9}, {1e9, -1e9}};
		for(int ss=idx; ss<=next_idx; ++ss) {
			const Real xBnd[2] = {myFish->rX[ss] - myFish->norX[ss]*myFish->width[ss],
					myFish->rX[ss] + myFish->norX[ss]*myFish->width[ss]};
			const Real yBnd[2] = {myFish->rY[ss] - myFish->norY[ss]*myFish->width[ss],
					myFish->rY[ss] + myFish->norY[ss]*myFish->width[ss]};
			const Real zBnd[2] = {-myFish->height[ss], +myFish->height[ss]};
			bbox[0][0] = std::min({bbox[0][0],xBnd[0],xBnd[1]});
			bbox[0][1] = std::max({bbox[0][1],xBnd[0],xBnd[1]});
			bbox[1][0] = std::min({bbox[1][0],yBnd[0],yBnd[1]});
			bbox[1][1] = std::max({bbox[1][1],yBnd[0],yBnd[1]});
			bbox[2][0] = std::min({bbox[2][0],zBnd[0],zBnd[1]});
			bbox[2][1] = std::max({bbox[2][1],zBnd[0],zBnd[1]});
		}
		vSegments[i].prepare(std::make_pair(idx, next_idx), bbox); //create a new segment
		vSegments[i].changeToComputationalFrame(position,quaternion);
	}

	// clear deformation velocities
	for(auto & entry : obstacleBlocks)
		delete entry.second;
	obstacleBlocks.clear();

	//MPI_Barrier(grid->getCartComm());
	t45 = std::chrono::high_resolution_clock::now();

	// 6. & 7.
	std::map<int, std::vector<VolumeSegment_OBB>> segmentsPerBlock;
	{
		for(int i=0;i<vInfo.size();++i) {
			const BlockInfo & info = vInfo[i];
			Real pStart[3], pEnd[3];
			info.pos(pStart, 0, 0, 0);
			info.pos(pEnd, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
			const Real safe_distance = 2.0*info.h_gridpoint; // two points on each side

			for(int s=0;s<Nsegments;++s)
				if(vSegments[s].isIntersectingWithAABB(pStart,pEnd,safe_distance))
					segmentsPerBlock[info.blockID].push_back(vSegments[s]);

			// allocate new blocks if necessary
			if(segmentsPerBlock.find(info.blockID) != segmentsPerBlock.end()) {
				assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
				obstacleBlocks[info.blockID] = new ObstacleBlock;
				obstacleBlocks[info.blockID]->clear();
			}
		}
	}
	//MPI_Barrier(grid->getCartComm());
	t67 = std::chrono::high_resolution_clock::now();

	assert(segmentsPerBlock.size() == obstacleBlocks.size());

	// 8.
	{
		#pragma omp parallel
		{
			PutFishOnBlocks putfish(myFish, position, quaternion);

			#pragma omp for schedule(dynamic)
			for(int i=0; i<vInfo.size(); i++) {
				BlockInfo info = vInfo[i];
				auto pos = segmentsPerBlock.find(info.blockID);
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;

				//tmpU will contain SDF: neg outside positive inside
				/*
					if(pos == segmentsPerBlock.end()) {
						for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
							for(int iy=0; iy<FluidBlock::sizeY; ++iy)
								for(int ix=0; ix<FluidBlock::sizeX; ++ix)
									b(ix,iy,iz).tmpU = -1.; //-1 here to avoid gremlins at blocks' boundaries
					} else {
				*/
				if(pos != segmentsPerBlock.end()) {
					for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
						for(int iy=0; iy<FluidBlock::sizeY; ++iy)
							for(int ix=0; ix<FluidBlock::sizeX; ++ix)
								b(ix,iy,iz).tmpU = 0.; //this will be accessed with plus equal

					assert(obstacleBlocks.find(info.blockID) != obstacleBlocks.end());
					ObstacleBlock* const defblock = obstacleBlocks.find(info.blockID)->second;
					putfish(info, b, defblock, pos->second);
				}
			}
		}
	}
	//MPI_Barrier(grid->getCartComm());
	t8 = std::chrono::high_resolution_clock::now();
/*
	{
		const auto len1   = std::chrono::duration<Real>(t1-t0).count();
		const auto len23  = std::chrono::duration<Real>(t23-t1).count();
		const auto len45  = std::chrono::duration<Real>(t45-t23).count();
		const auto len67  = std::chrono::duration<Real>(t67-t45).count();
		const auto len8   = std::chrono::duration<Real>(t8-t67).count();
		printf("Creation times %d %d: %g %g %g %g %g \n",rank, obstacleID, len1, len23, len45, len67, len8);
		fflush(0);
	}
*/
	ptrUinf_copy = Uinf;
}

void IF3D_FishOperator::finalize(const int step_id,const Real time, const Real dt, const Real *Uinf)
{
	// STRATEGY
	// 9. create the Chi out of the SDF. In same sweep, compute the actual CoM
	// 10. compute all shit: linear momentum, angular momentum etc.
	// 11. correct deformation velocity to nullify momenta for the final discrete representation
	std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1, t23, t45, t67, t8, t910, t11;

	t8 = std::chrono::high_resolution_clock::now();
	// 9. & 10.
	{
		const int nthreads = omp_get_max_threads();
		vector<surfaceBlocks> dataPerThread(nthreads);
		vector<array<Real,4>> momenta(nthreads);
		for(int i=0; i<nthreads; i++) for(int j=0; j<4; j++) momenta[i][j]=0;

		vector<PutFishOnBlocks_Finalize*> finalize;
		for(int i=0; i<nthreads; i++) {
			PutFishOnBlocks_Finalize* tmp = new
        PutFishOnBlocks_Finalize(&obstacleBlocks,&dataPerThread[i],&momenta[i]);
			finalize.push_back(tmp);
		}
		compute(finalize);
		for(int i=0; i<nthreads; i++) delete finalize[i];

		double sumX[4] = {0,0,0,0};
		double totX[4] = {0,0,0,0};
		for(int i=0; i<nthreads; i++) {
			sumX[0] += momenta[i][0];
			sumX[1] += momenta[i][1];
			sumX[2] += momenta[i][2];
			sumX[3] += momenta[i][3];
		}

		MPI_Allreduce(sumX, totX, 4, MPI::DOUBLE, MPI::SUM, grid->getCartComm());
		surfData.finalizeOnGrid(dataPerThread);

		assert(totX[0]>std::numeric_limits<double>::epsilon());
		CoM_interpolated[0]=totX[1]/totX[0];
		CoM_interpolated[1]=totX[2]/totX[0];
		CoM_interpolated[2]=totX[3]/totX[0];
	}
	//MPI_Barrier(grid->getCartComm());
	t910 = std::chrono::high_resolution_clock::now();

	#ifdef __useSkin_
	myFish->surfaceToComputationalFrame(_2Dangle,CoM_interpolated);
	#endif
	// 11.
	_makeDefVelocitiesMomentumFree(CoM_interpolated);
	//MPI_Barrier(grid->getCartComm());
	t11 = std::chrono::high_resolution_clock::now();
/*
	{
		const auto len910 = std::chrono::duration<Real>(t910-t8).count();
		const auto len11  = std::chrono::duration<Real>(t11-t910).count();
		printf("Finalization times %d %d: %g %g\n",rank, obstacleID, len910, len11);
		fflush(0);
	}
*/
}

void IF3D_FishOperator::interpolateOnSkin(const Real time, const int stepID, bool _dumpWake)
{
	#ifdef __useSkin_
  assert(quaternion[1] == 0 && quaternion[2] == 0);
	sr.updateStepId(stepID+obstacleID);
	myFish->computeSkinNormals(_2Dangle, CoM_interpolated);
  sr.nearestGridPoints(&surfData, myFish->upperSkin->Npoints,
                        myFish->upperSkin->xSurf, myFish->upperSkin->ySurf,
                        myFish->lowerSkin->xSurf, myFish->lowerSkin->ySurf,
                        myFish->upperSkin->normXSurf, myFish->upperSkin->normYSurf,
                        myFish->lowerSkin->normXSurf, myFish->lowerSkin->normYSurf,
                        position[2], vInfo[0].h_gridpoint, grid->getCartComm());

		#ifndef __RL_TRAINING
			if(rank==0) sr.print(obstacleID, stepID, time);
  	#endif

	if(_dumpWake && ptrUinf_copy not_eq nullptr)
			dumpWake(stepID, time, ptrUinf_copy);
	if(ptrUinf_copy == nullptr && !rank) printf("(null backup of uinf)\n");
  #endif
}

void IF3D_FishOperator::update(const int stepID, const Real t, const Real dt, const Real *Uinf)
{
	// synchronize internal time
	sim_time = t + dt;
	sim_dt = dt;
	// update position and angles
	IF3D_ObstacleOperator::update(stepID,t, dt, Uinf);
	// negative: we subtracted this angvel
	theta_internal -= sim_dt*angvel_internal;
	angvel_integral[0] += dt*angVel[0];
	angvel_integral[1] += dt*angVel[1];
	angvel_integral[2] += dt*angVel[2];
}

void IF3D_FishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_ObstacleOperator::_parseArguments(parser);
	parser.set_strict_mode();
	parser.unset_strict_mode();
	Tperiod = parser("-T").asDouble(1.0);
	nActions = parser("-nActions").asInt(0);
	GoalDX = parser("-GoalDX").asDouble(2.0);
	phaseShift = parser("-phi").asDouble(0.0);
	Tstartlearn = parser("-Tstartlearn").asDouble(1e6);
	bCorrectTrajectory = parser("-Correct").asBool(false);
	bInteractive = parser("-Active").asBool(false);
	NpLatLine = parser("-NpLatLine").asInt(__NpLatLine);
	followX = parser("-followX").asDouble(-1);
	followY = parser("-followY").asDouble(-1);
	if(NpLatLine != __NpLatLine) {
		printf("Mismatch in __NpLatLine, check settings\n");
		fflush(0); abort();
	}
	sr.set_NpLatLine(NpLatLine);
	//i want to reset time-averaged quantities before first actual comm
	sr.t_next_comm = Tstartlearn - Tperiod/2.;
	sr.bForgiving = parser("-easyFailBox").asBool(false);
	sr.GoalDX = GoalDX;
  sr.thExp = _2Dangle;
	#ifdef __useSkin_
	bHasSkin = true;
	#endif
	/*
		randomStart = parser("-randomStart").asBool(false);
		if (randomStart) {
			printf("Random start\n");
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<Real> dis(-1.,1.);
			position[0] += .5*length*dis(gen);
			position[1] += .1*length*dis(gen);
		}
	 */
}

void IF3D_FishOperator::getCenterOfMass(Real CM[3]) const
{
	// return computation CoM, not the one were advecting
	CM[0]=CoM_interpolated[0];
	CM[1]=CoM_interpolated[1];
	CM[2]=CoM_interpolated[2];
}

void IF3D_FishOperator::getSkinsAndPOV(Real& x, Real& y, Real& th,
  Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts)
{
	assert(quaternion[1] == 0 && quaternion[2] == 0);
	x  = position[0];
	y  = position[1];
	th  = _2Dangle;
	pXL = myFish->lowerSkin->xSurf;
	pYL = myFish->lowerSkin->ySurf;
	pXU = myFish->upperSkin->xSurf;
	pYU = myFish->upperSkin->ySurf;
	Npts = myFish->lowerSkin->Npoints;
}

void IF3D_FishOperator::save(const int stepID, const Real t, string filename)
{
    //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename + ".txt");

    savestream<<t<<"\t"<<sim_dt<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
    savestream<<_2Dangle;
    savestream.close();

}

void IF3D_FishOperator::restart(const Real t, string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream >> _2Dangle;
    restartstream.close();

	std::cout<<"RESTARTED FISH: "<<std::endl;
	std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
	std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
	std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
	std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
	std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
	std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
	std::cout<<"2D angle: \t"<<_2Dangle<<std::endl;
}

/*
void IF3D_CarlingFishOperator::writeToFile(const int step_id, const Real t, std::string filename)
{
    std::string fname = (filename==std::string()) ? "fish" : filename;

    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);

    savestream.open(fname+"_interpolated.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "CoM[2]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "linMom[2]\t";
        savestream << "angMom[0]\t";
        savestream << "angMom[1]\t";
        savestream << "angMom[2]" << std::endl;
    }

    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << object_ongrid.volume << "\t";
    savestream << CoM_interpolated[0] << "\t";
    savestream << CoM_interpolated[1] << "\t";
    savestream << CoM_interpolated[2] << "\t";
    savestream << object_ongrid.linearMomentum[0] << "\t";
    savestream << object_ongrid.linearMomentum[1] << "\t";
    savestream << object_ongrid.linearMomentum[2] << "\t";
    savestream << object_ongrid.angularMomentum[0] << "\t";
    savestream << object_ongrid.angularMomentum[1] << "\t";
    savestream << object_ongrid.angularMomentum[2] << "\t";
    savestream << object_ongrid.J[2] << std::endl;
    savestream.close();

    savestream.open(fname+"_internal.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "angMom\t";
        savestream << "theta\t";
        savestream << "angvel" << std::endl;
    }

    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << volume_internal << "\t";
    savestream << CoM_internal[0] << "\t";
    savestream << CoM_internal[1] << "\t";
    savestream << vCoM_internal[0]*volume_internal << "\t";
    savestream << vCoM_internal[1]*volume_internal << "\t";
    savestream << angvel_internal*J_internal << "\t";
    savestream << theta_internal << "\t";
    savestream << angvel_internal << std::endl;
    savestream.close();

    savestream.open(fname+"_computation.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "pos[0]\t";
        savestream << "pos[1]\t";
        savestream << "pos[2]\t";
        savestream << "quat[0]\t";
        savestream << "quat[1]\t";
        savestream << "quat[2]\t";
        savestream << "quat[3]\t";
        savestream << "transVel[0]\t";
        savestream << "transVel[1]\t";
        savestream << "transVel[2]\t";
        savestream << "angVel[0]\t";
        savestream << "angVel[1]\t";
        savestream << "angVel[2]" << std::endl;
    }

    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << position[0] << "\t";
    savestream << position[1] << "\t";
    savestream << position[2] << "\t";
    savestream << quaternion[0] << "\t";
    savestream << quaternion[1] << "\t";
    savestream << quaternion[2] << "\t";
    savestream << quaternion[3] << "\t";
    savestream << transVel[0] << "\t";
    savestream << transVel[1] << "\t";
    savestream << transVel[2] << "\t";
    savestream << angVel[0] << "\t";
    savestream << angVel[1] << "\t";
    savestream << angVel[2] << std:: endl;
    savestream.close();
}*/
