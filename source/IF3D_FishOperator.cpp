//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_FishOperator.h"
#include "IF3D_FishLibrary.h"

Fish::FishMidlineData::FishMidlineData(const int Nm, const Real len, const Real Tp, const Real phase, const Real dx_ext):
Nm(Nm),length(len),Tperiod(Tp),phaseShift(phase),rS(_alloc(Nm)),rX(_alloc(Nm)),rY(_alloc(Nm)),
vX(_alloc(Nm)),vY(_alloc(Nm)),norX(_alloc(Nm)),norY(_alloc(Nm)),vNorX(_alloc(Nm)),vNorY(_alloc(Nm)),
width(_alloc(Nm)),height(_alloc(Nm)),iFishStart(4*NPPEXT),iFishEnd(Nm-1-4*NPPEXT)
{
	// extension_info contains number of extension points and extension dx
	const int Nextension = 4*NPPEXT; // up to 3dx on each side (to get proper interpolation up to 2dx)
	const int Next = Nextension; // number of points per extension
	const int Nint = Nm -2*Next; // number of interior points

	// extension head
	for(int i=0;i<Next;++i)
		rS[i] = 0.0 - (Next- i) * dx_ext;
	// interior points
	for(int i=0;i<Nint;++i)
		rS[i+Next] = length * 0.5 * (1.0 - std::cos(i * M_PI/((Real)Nint-1))); // cosine: more points near head and tail
	// rS[i] = i*length/((Real)Nint-1); // linear: equally distributed points
	// extension tail
	for(int i=0;i<Next;++i)
		rS[i+Nint+Next] = length + (i + 1)*dx_ext;
	_computeWidthsHeights();
}

Fish::FishMidlineData::~FishMidlineData()
{
    _dealloc(rS);
    _dealloc(rX);
    _dealloc(rY);
    _dealloc(vX);
    _dealloc(vY);
    _dealloc(norX);
    _dealloc(norY);
    _dealloc(vNorX);
    _dealloc(vNorY);
    _dealloc(height);
    _dealloc(width);
}

void Fish::FishMidlineData::_prepareRotation2D(const Real angle)
{
	Rmatrix2D[0][0] = Rmatrix2D[1][1] = std::cos(angle);
	Rmatrix2D[0][1] = -std::sin(angle);
	Rmatrix2D[1][0] = -Rmatrix2D[0][1];
}

void Fish::FishMidlineData::_computeWidthsHeights()
{
	for(int i=0;i<Nm;++i) {
		width[i]  = Fish::width(rS[i],length);
		height[i] = Fish::height(rS[i],length);
	}
}

void Fish::FishMidlineData::_computeMidlineNormals()
{
#pragma omp parallel for
	for(int i=0; i<Nm-1; i++) {
		const Real ds = rS[i+1]-rS[i];
		const Real tX = rX[i+1]-rX[i];
		const Real tY = rY[i+1]-rY[i];
		const Real tVX = vX[i+1]-vX[i];
		const Real tVY = vY[i+1]-vY[i];
		norX[i] = -tY/ds;
		norY[i] =  tX/ds;
		vNorX[i] = -tVY/ds;
		vNorY[i] =  tVX/ds;
	}
	norX[Nm-1] = norX[Nm-2];
	norY[Nm-1] = norY[Nm-2];
	vNorX[Nm-1] = vNorX[Nm-2];
	vNorY[Nm-1] = vNorY[Nm-2];
}

Real Fish::FishMidlineData::integrateLinearMomentum(Real CoM[2], Real vCoM[2])
{   // already worked out the integrals for r, theta on paper
	// remaining integral done with composite trapezoidal rule
	// minimize rhs evaluations --> do first and last point separately
	Real _vol(0), _cmx(0), _cmy(0), _lmx(0), _lmy(0);
#pragma omp parallel for reduction(+:_vol,_cmx,_cmy,_lmx,_lmy)
	for(int i=0;i<Nm;++i) {
		const Real ds = (i==0) ? rS[1]-rS[0] :
				((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
		const Real fac1 = _integrationFac1(i);
		const Real fac2 = _integrationFac2(i);
		_vol += 0.5*fac1*ds;
		_cmx += 0.5*(rX[i]*fac1 + norX[i]*fac2)*ds;
		_cmy += 0.5*(rY[i]*fac1 + norY[i]*fac2)*ds;
		_lmx += 0.5*(vX[i]*fac1 + vNorX[i]*fac2)*ds;
		_lmy += 0.5*(vY[i]*fac1 + vNorY[i]*fac2)*ds;
	}

	vol=_vol*M_PI;
	CoM[0]=_cmx*M_PI;
	CoM[1]=_cmy*M_PI;
	linMom[0]=_lmx*M_PI;
	linMom[1]=_lmy*M_PI;

	assert(vol> std::numeric_limits<Real>::epsilon());
	const Real ivol = 1.0/vol;

	CoM[0]*=ivol;
	CoM[1]*=ivol;
	vCoM[0]=linMom[0]*ivol;
	vCoM[1]=linMom[1]*ivol;
	//printf("%f %f %f %f %f\n",CoM[0],CoM[1],vCoM[0],vCoM[1], vol);
	return vol;
}

void Fish::FishMidlineData::integrateAngularMomentum(Real & angVel)
{
	// assume we have already translated CoM and vCoM to nullify linear momentum

	// already worked out the integrals for r, theta on paper
	// remaining integral done with composite trapezoidal rule
	// minimize rhs evaluations --> do first and last point separately
	Real _J(0), _am(0);

#pragma omp parallel for reduction(+:_J,_am)
	for(int i=0;i<Nm;++i) {
		const Real ds = (i==0) ? rS[1]-rS[0] :
				((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
		const Real fac1 = _integrationFac1(i);
		const Real fac2 = _integrationFac2(i);
		const Real fac3 = _integrationFac3(i);
		double tmp_J, tmp_M;
		tmp_M  = (rX[i]*vY[i] - rY[i]*vX[i])*fac1;
		tmp_M += (rX[i]*vNorY[i] - rY[i]*vNorX[i] + vY[i]*norX[i] - vX[i]*norY[i])*fac2;
		tmp_M += (norX[i]*vNorY[i] - norY[i]*vNorX[i])*fac3;
		_am += 0.5*tmp_M*ds;
		tmp_J  = (rX[i]*rX[i] + rY[i]*rY[i])*fac1;
		tmp_J += 2.0*(rX[i]*norX[i] + rY[i]*norY[i])*fac2;
		//tmpSum += (norX[idx]*norX[idx]+norY[idx]*norY[idx])*fac3;
		tmp_J += fac3;
		_J += 0.5*tmp_J*ds;
	}

	J=_J*M_PI;
	angMom=_am*M_PI;
	assert(J>std::numeric_limits<Real>::epsilon());
	angVel = angMom/J;
}

void Fish::FishMidlineData::changeToCoMFrameLinear(const Real CoM_internal[2], const Real vCoM_internal[2])
{
	for(int i=0;i<Nm;++i) {
		rX[i]-=CoM_internal[0];
		rY[i]-=CoM_internal[1];
		vX[i]-=vCoM_internal[0];
		vY[i]-=vCoM_internal[1];
	}
}

void Fish::FishMidlineData::changeToCoMFrameAngular(const Real theta_internal, const Real angvel_internal)
{
	_prepareRotation2D(theta_internal);
#pragma omp parallel for
	for(int i=0;i<Nm;++i) {
		_rotate2D(rX[i],rY[i]);
		_rotate2D(vX[i],vY[i]);
		vX[i] += angvel_internal*rY[i];
		vY[i] -= angvel_internal*rX[i];
	}
	_computeMidlineNormals();
}

IF3D_FishOperator::IF3D_FishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_ObstacleOperator(grid, parser), theta_internal(0.0), angvel_internal(0.0), sim_time(0.0),
  sim_dt(0.0), adjTh(adjTh), myFish(nullptr), angvel_integral{0.,0.,0.}
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

	const int Nsegments = NPPSEG;
	const int Nextension = 4*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
	const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
	const Real target_ds = vInfo[0].h_gridpoint/TGTPPB;
	const Real target_Nm = length/target_ds;
	const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG) + 1;
	assert((Nm-1)%Nsegments==0);
	if (bCorrectTrajectory) {
		Real velx_tot = Uinf[0] - transVel[0];
		Real vely_tot = Uinf[1] - transVel[1];
		Real AngDiff  = std::atan2(vely_tot,velx_tot);
		adjTh = (1.-dt) * adjTh + dt * AngDiff;
		const Real B = (AngDiff*angVel[2]>0) ? 0.25/M_PI : 0;
		const Real PID = .5*adjTh +B*AngDiff*fabs(angVel[2]);
		myFish->_correctTrajectory(PID, time, dt);
	}
	// 1.
	myFish->computeMidline(time);

	// 2. & 3.
	volume_internal = myFish->integrateLinearMomentum(CoM_internal, vCoM_internal);
	assert(volume_internal > std::numeric_limits<Real>::epsilon());
	myFish->changeToCoMFrameLinear(CoM_internal, vCoM_internal);

	angvel_internal_prev = angvel_internal;
	myFish->integrateAngularMomentum(angvel_internal);
	J_internal = myFish->J;
	// update theta now with new angvel info
	//theta_internal -= 0.5*sim_dt*(angvel_internal+angvel_internal_prev);//negative: we subtracted this angvel
	myFish->changeToCoMFrameAngular(theta_internal, angvel_internal);

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

	// 4.
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

	//assert(not segmentsPerBlock.empty()); //killed this assert: distributed fish
	assert(segmentsPerBlock.size() == obstacleBlocks.size());

	// 8.
	{
#pragma omp parallel
		{
			PutFishOnBlocks putfish(myFish, position, quaternion);

#pragma omp for schedule(static)
			for(int i=0; i<vInfo.size(); i++) {
				BlockInfo info = vInfo[i];
				auto pos = segmentsPerBlock.find(info.blockID);
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;

				//tmpU will contain SDF: neg outside positive inside
				if(pos == segmentsPerBlock.end()) {
					for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
						for(int iy=0; iy<FluidBlock::sizeY; ++iy)
							for(int ix=0; ix<FluidBlock::sizeX; ++ix)
								b(ix,iy,iz).tmpU = -1.; //-1 here to avoid gremlins at blocks' boundaries
				} else {
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

	// 9. & 10. & 11.
	{
		const int nthreads = omp_get_max_threads();
		vector<surfaceBlocks> dataPerThread(nthreads);
		vector<array<Real,4>> momenta(nthreads);
		vector<PutFishOnBlocks_Finalize*> finalize;
		for(int i = 0; i < nthreads; ++i) {
			PutFishOnBlocks_Finalize* tmp = new PutFishOnBlocks_Finalize(&obstacleBlocks,&dataPerThread[i],&momenta[i]);
			finalize.push_back(tmp);
		}
		compute(finalize);

		double sumX[4] = {0,0,0,0};
		double totX[4] = {0,0,0,0};
		for(int i=0; i<nthreads; i++) {
			sumX[0] += momenta[i][0];
			sumX[1] += momenta[i][1];
			sumX[2] += momenta[i][2];
			sumX[3] += momenta[i][3];
		}

		MPI::COMM_WORLD.Allreduce(sumX, totX, 4, MPI::DOUBLE, MPI::SUM);

		surfData.finalizeOnGrid(dataPerThread);

		assert(totX[0]>std::numeric_limits<double>::epsilon());
		CoM_interpolated[0]=totX[1]/totX[0];
		CoM_interpolated[1]=totX[2]/totX[0];
		CoM_interpolated[2]=totX[3]/totX[0];

		_makeDefVelocitiesMomentumFree(CoM_interpolated);
	}
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
	Tperiod = parser("-T").asDouble();
	parser.unset_strict_mode();
	nActions = parser("-nActions").asInt(0);
	GoalDX = parser("-GoalDX").asDouble(0.0);
	phaseShift = parser("-phi").asDouble(0.0);
	Tstartlearn = parser("-Tstartlearn").asDouble(1e6);
	bCorrectTrajectory = parser("-Correct").asBool(false);
    randomStart = parser("-randomStart").asBool(false);
    if (randomStart) {
    	printf("Random start\n");
    	std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<Real> dis(-1.,1.);
    	position[0] += .5*length*dis(gen);
    	position[1] += .1*length*dis(gen);
    }
    /*
    //TODO state and reward:
    sr->updateInstant(position[0], position[1], angle, 0., 0., 0.);
    sr->t_next_comm = Tstartlearn - 1/2.; //i want to reset time-averages before first actual comm
    bool bForgiving = parser("-easyFailBox").asBool(false);
    sr->bForgiving = bForgiving;
    sr->GoalDX = GoalDX;
    sr->thExp = angle;
    */
}

/*
void IF3D_StefanLearnTurnOperator::execute(Communicator * comm, const int iAgent, const Real time)
{
    if (time < Tstartlearn) {
        sr->resetAverage();
        sr->t_next_comm = Tstartlearn;

        //TMP: first rnd action is going to be taken after a while
        if (randomActions) sr->t_next_comm = Tstartlearn+6.;

        return;
    }

    if (not bInteractive) {
        if (not randomActions) { sr->t_next_comm=1e6; return; }
        //we might decide to pause turning, then just pick a pause counter
        if (nPauseActions-- > 0) {
            vector<Real> raction(1,0.);
            myFish->execute(time, sr->t_next_comm, raction);
            sr->t_next_comm += .5*myFish->Tperiod;
            printf("pausing from turning at time %f, will still pause %d turns.
            		Next turn at time %f\n",time,nPauseActions,sr->t_next_comm);
            return;
        }
        vector<Real> raction(1,signLastTurn);
        myFish->execute(time, sr->t_next_comm, raction);
        sr->t_next_comm += .5*myFish->Tperiod;
        printf("turning at time %f with modifier %f, turn counter is %d.
        		Next turn at time %f\n",time,signLastTurn,nTurnActions,sr->t_next_comm);

        if (++nTurnActions >= 4) {
            nPauseActions = 10.;
            nTurnActions = 0;
            signLastTurn *= -1.;
        }

    } else if (useLoadedActions) {

        vector<Real> actions(nActions);
        if (loadedActions.size()>1) {
            actions = loadedActions.back();
            loadedActions.pop_back();
        } //else zero actions
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }
        sr->resetAverage();

    } else {

        const Real relT= fmod(time,1.); //1 is Tperiod
#ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*NpLatLine : 25+  8*NpLatLine;
#else
        const int nStates = (nActions==1) ? 20+10*NpLatLine : 25+ 10*NpLatLine;
#endif
        vector<Real> state(nStates), actions(nActions);

        int k(0);
        state[k++] = sr->Xrel - GoalDX;
        state[k++] = sr->Yrel;
        state[k++] = sr->RelAng;
        state[k++] = relT;
        state[k++] = new_curv;
        state[k++] = old_curv;

        if(nActions==2) { //this is for backwards compatibility
            state[k++] = new_Tp;
                        //2.*M_PI*((time-time0)/l_Tp +timeshift -rS[i]/length) + M_PI*phaseShift
            Real Fshift = 2.*((-myFish->time0)/myFish->l_Tp +myFish->timeshift)+myFish->phaseShift;
            Fshift = fmod(Fshift,2.0);
            state[k++] = (Fshift<0) ? 2.+Fshift : Fshift;
            state[k++] = sr->VX;
            state[k++] = sr->VY;
            state[k++] = sr->AV;
        }

        state[k++] = sr->Dist;
        state[k++] = sr->Quad;
        state[k++] = sr->VxAvg;
        state[k++] = sr->VyAvg;
        state[k++] = sr->AvAvg;
        state[k++] = sr->Pout;
        state[k++] = sr->defPower;
        state[k++] = sr->EffPDef;
        state[k++] = sr->PoutBnd;
        state[k++] = sr->defPowerBnd;
        state[k++] = sr->EffPDefBnd;
        state[k++] = sr->Pthrust;
        state[k++] = sr->Pdrag;
        state[k++] = sr->ToD;

        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVBelow[j];
        #ifndef _NOVISION_
        for (int j=0; j<2*NpLatLine; j++) state[k++] = sr->raySight[j];
        #endif
        const Real reward = (sr->info==2) ? -10 : sr->EffPDefBnd;
        comm->sendState(iAgent-1, sr->info, state, reward); //TODO
        fflush(0);
        if (sr->info==2) return;

        sr->info = 0;

        comm->recvAction(actions);
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }

        #ifndef TRAINING
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        if(nActions==2)
            filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
        #endif //TRAINING

        sr->resetAverage();
    }
}
*/

void IF3D_FishOperator::getCenterOfMass(Real CM[3]) const
{
	// return computation CoM, not the one were advecting
	CM[0]=CoM_interpolated[0];
	CM[1]=CoM_interpolated[1];
	CM[2]=CoM_interpolated[2];
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
    restartstream.close();
    
	std::cout<<"RESTARTED FISH: "<<std::endl;
	std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
	std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
	std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
	std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
	std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
	std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
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
