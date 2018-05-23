//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_FishOperator.h"
#include "IF3D_FishLibrary.h"
#include <chrono>
IF3D_FishOperator::IF3D_FishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) : IF3D_ObstacleOperator(g, p, u)
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

aryVolSeg IF3D_FishOperator::prepare_vSegments(const int Nsegments, const int Nm)
{
  aryVolSeg vSegments(Nsegments);
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
      const Real maxX=std::max(xBnd[0],xBnd[1]), minX=std::min(xBnd[0],xBnd[1]);
      const Real maxY=std::max(yBnd[0],yBnd[1]), minY=std::min(yBnd[0],yBnd[1]);
      const Real maxZ=std::max(zBnd[0],zBnd[1]), minZ=std::min(zBnd[0],zBnd[1]);
      bbox[0][0] = std::min(bbox[0][0], minX);
      bbox[0][1] = std::max(bbox[0][1], maxX);
      bbox[1][0] = std::min(bbox[1][0], minY);
      bbox[1][1] = std::max(bbox[1][1], maxY);
      bbox[2][0] = std::min(bbox[2][0], minZ);
      bbox[2][1] = std::max(bbox[2][1], maxZ);
    }
    vSegments[i].prepare(std::make_pair(idx, next_idx), bbox); //create a new segment
    vSegments[i].changeToComputationalFrame(position,quaternion);
  }
  return vSegments;
}

mapBlock2Segs IF3D_FishOperator::prepare_segPerBlock(const aryVolSeg&vSegments)
{
  mapBlock2Segs segmentsPerBlock;

  // clear deformation velocities
  for(auto & entry : obstacleBlocks) delete entry.second;
  obstacleBlocks.clear();

  for(int i=0; i<vInfo.size(); ++i) {
    const BlockInfo & info = vInfo[i];
    Real pStart[3], pEnd[3];
    info.pos(pStart, 0, 0, 0);
    info.pos(pEnd, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    const Real safe_distance = 2.0*info.h_gridpoint; // two points on each side
    //const Real safe_distance = info.h_gridpoint; // one point on each side for Towers

    for(int s=0; s<vSegments.size(); ++s)
      if(vSegments[s].isIntersectingWithAABB(pStart,pEnd,safe_distance))
        segmentsPerBlock[info.blockID].push_back(vSegments[s]);

    // allocate new blocks if necessary
    if(segmentsPerBlock.find(info.blockID) != segmentsPerBlock.end()) {
      assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
      ObstacleBlock * const block = new ObstacleBlock();
      assert(block not_eq nullptr);
      obstacleBlocks[info.blockID] = block;
      block->clear();
    }
  }
  return segmentsPerBlock;
}

void IF3D_FishOperator::apply_pid_corrections(const double time, const double dt, const Real *Uinf)
{
  if (bCorrectTrajectory && time>0.312)
  {
    assert(followX < 0 && followY < 0);
    const double velx_tot = Uinf[0] - transVel[0];
    const double vely_tot = Uinf[1] - transVel[1];
    const double AngDiff  = std::atan2(vely_tot,velx_tot);
    adjTh = (1.-dt) * adjTh + dt * AngDiff;
    const double INST = (AngDiff*angVel[2]>0) ? AngDiff*std::fabs(angVel[2]) : 0;
    const double PID = 0.1*adjTh + 0.01*INST;
    myFish->_correctTrajectory(PID, 0,time, dt);
  }

  if (followX > 0 && followY > 0) //then i control the position
  {
    assert(not bCorrectTrajectory);
    const double velx_tot = Uinf[0] - transVel[0];
    const double vely_tot = Uinf[1] - transVel[1];
    const double AngDiff  = _2Dangle;//std::atan2(vely_tot,velx_tot);

    // Control posDiffs
    const double xDiff = (position[0] - followX)/length;
    const double yDiff = (position[1] - followY)/length;
    const double absDY = std::fabs(yDiff);
    const double velAbsDY = yDiff>0 ? transVel[1]/length : -transVel[1]/length;
    const double velDAvg = AngDiff-adjTh + dt*angVel[2];

    adjTh = (1.-dt) * adjTh + dt * AngDiff;
    adjDy = (1.-dt) * adjDy + dt * yDiff;

    //If angle is positive: positive curvature only if Dy<0 (must go up)
    //If angle is negative: negative curvature only if Dy>0 (must go down)
    //const Real INST = (AngDiff*angVel[2]>0 && yDiff*AngDiff<0) ? AngDiff*std::fabs(yDiff*angVel[2]) : 0;
    const double PROP = (adjTh  *yDiff<0) ?   adjTh*absDY : 0;
    const double INST = (AngDiff*yDiff<0) ? AngDiff*absDY : 0;

    //zero also the derivatives when appropriate
    const double f1=(PROP ? 20 : 0), f2=(INST ? 50 : 0), f3=1;

    // Linearly increase (or decrease) amplitude to 1.2X (decrease to 0.8X)
    //(experiments observed 1.2X increase in amplitude when swimming faster)
    //if fish falls back 1 body length. Beyond that, will still increase but dunno if will work
    const double ampFac = f3*xDiff + 1.0;
    const double ampVel = f3*transVel[0]/length;

    const double curv1fac = f1*PROP;
    const double curv1vel = f1*(velAbsDY*adjTh   + absDY*velDAvg);
    const double curv2fac = f2*INST;
    const double curv2vel = f2*(velAbsDY*AngDiff + absDY*angVel[2]);
                //const Real vPID = velAbsDY*(f1*adjTh + f2*AngDiff) + absDY*(f1*velDAvg+f2*angVel[2]);
                //const Real PID = f1*PROP + f2*INST;
    if(!rank) printf("%f\t f1: %f %f\t f2: %f %f\t f3: %f %f\n", time,
      curv1fac, curv1vel, curv2fac, curv2vel, ampFac, ampVel);
    myFish->_correctTrajectory(curv1fac+curv2fac, curv1vel+curv2vel, time, dt);
    myFish->_correctAmplitude(ampFac, ampVel, time, dt);
  }
}

void IF3D_FishOperator::integrateMidline()
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

  #ifndef NDEBUG
  {
    double dummy_CoM_internal[2], dummy_vCoM_internal[2], dummy_angvel_internal;
    // check that things are zero
    const double volume_internal_check = myFish->integrateLinearMomentum(dummy_CoM_internal,dummy_vCoM_internal);
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
}

void IF3D_FishOperator::writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock)
{
  #pragma omp parallel
  {
    PutFishOnBlocks putfish(myFish, position, quaternion);

    #pragma omp for schedule(dynamic)
    for(int i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      const auto pos = segmentsPerBlock.find(info.blockID);
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;

      if(pos != segmentsPerBlock.end()) {
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          b(ix,iy,iz).tmpU = 0.; //this will be accessed with plus equal

        assert(obstacleBlocks.find(info.blockID) != obstacleBlocks.end());
        ObstacleBlock*const block = obstacleBlocks.find(info.blockID)->second;
        putfish(info, b, block, pos->second);
      }
    }
  }
}

void IF3D_FishOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{
  // STRATEGY
  // we need some things already
  // - the internal angle at the previous timestep, obtained from integrating the actual def velocities
  // (not the imposed deformation velocies, because they dont have zero ang mom)
  // - the internal angular velocity at previous timestep

  // 1. create midline
  // 2. integrate to find CoM, angular velocity, etc
  // 3. shift midline to CoM frame: zero internal linear momentum and angular momentum

  // 4. split the fish into segments (according to s)
  // 5. rotate the segments to computational frame (comp CoM and angle)
  // 6. for each Block in the domain, find those segments that intersect it
  // 7. for each of those blocks, allocate an ObstacleBlock

  // 8. put the 3D shape on the grid: SDF-P2M for sdf, normal P2M for udef

  const int Nextension = NEXTDX*NPPEXT; // up to 3dx on each side (to get proper interpolation up to 2dx)
  const double dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
  const double target_ds = vInfo[0].h_gridpoint/TGTPPB;
  const double target_Nm = length/target_ds;
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

  apply_pid_corrections(time, dt, Uinf);

  // 1.
  myFish->computeMidline(time);
  #ifdef __useSkin_
  myFish->computeSurface();
  #endif

  // 2. & 3.
  integrateMidline();

  //CAREFUL: this func assumes everything is already centered around CM to start with, which is true (see steps 2. & 3. ...) for rX, rY: they are zero at CM, negative before and + after

  // To output the raw midlines
  #if 0
  // Check if we must take a dump
  double currentDumpFactor = time/0.0001;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if((std::floor(currentDumpFactor) > previousDumpFactor) & rank==0 )
  {
    //dumpNow = true;
    char buf[500];
    sprintf(buf, "midline_%07d.txt", step_id);
    FILE * f = fopen(buf,"w");
    fprintf(f, "s x y xGlob yGlob uBody vBody\n");
    for (int i=0; i<myFish->Nm; i++){
     double temp[3] = {myFish->rX[i], myFish->rY[i], 0.0};
     dummy.changeToComputationalFrame(temp);
     Real udef[3] = {myFish->vX[i], myFish->vY[i], 0.0};
     dummy.changeVelocityToComputationalFrame(udef);
     fprintf(f, "%g %g %g %g %g %g %g\n",
     myFish->rS[i],myFish->rX[i],myFish->rY[i],temp[0],temp[1],udef[0],udef[1]);
    }
    printf("Dumped midline\n");
    fclose(f);
  }
  previousDumpFactor = floor(currentDumpFactor);
  #endif

  // 4. & 5.
  const auto vSegments = prepare_vSegments(Nsegments, Nm);

  // 6. & 7.
  const auto segmentsPerBlock = prepare_segPerBlock(vSegments);
  assert(segmentsPerBlock.size() == obstacleBlocks.size());

  // 8.
  writeSDFOnBlocks(segmentsPerBlock);
}

void IF3D_FishOperator::computeChi(const int step_id, const double time, const double dt, const Real*Uinf, int& mpi_status)
{
  // 9. create the Chi out of the SDF. In same sweep, compute the actual CoM
  if(obstacleBlocks.size() || mpi_status==2)
  {
    if (mpi_status==1)
    {
      printf("Each rank can call computeChi only once. Fish got too close\n");
      fflush(0);
      MPI_Abort(grid->getCartComm(), 0);
    }

    const int nthreads = omp_get_max_threads();

    vector<PutFishOnBlocks_Finalize*> finalize(nthreads, nullptr);
    for(int i=0; i<nthreads; i++) finalize[i] = new PutFishOnBlocks_Finalize();

    compute<PutFishOnBlocks_Finalize, VOLUME>(finalize);

    for(int i=0; i<nthreads; i++) delete finalize[i];
    mpi_status = 1;
  }
}

void IF3D_FishOperator::finalize(const int step_id,const double time, const double dt, const Real *Uinf)
{
  // 10. compute all shit: linear momentum, angular momentum etc.
  // 11. correct deformation velocity to nullify momenta for the final discrete representation

  // 10.
  vector<double> com(4, 0);
  for (auto & block : obstacleBlocks) {
    com[0] += block.second->mass;
    com[1] += block.second->CoM_x;
    com[2] += block.second->CoM_y;
    com[3] += block.second->CoM_z;
  }

  MPI_Allreduce(MPI_IN_PLACE, com.data(), 4, MPI::DOUBLE, MPI::SUM, grid->getCartComm());

  assert(com[0]>std::numeric_limits<Real>::epsilon());
  CoM_interpolated[0]=com[1]/com[0];
  CoM_interpolated[1]=com[2]/com[0];
  CoM_interpolated[2]=com[3]/com[0];

  #ifdef __useSkin_
  myFish->surfaceToComputationalFrame(_2Dangle,CoM_interpolated);
  #endif

  // 11.
  _makeDefVelocitiesMomentumFree(CoM_interpolated);
  for(auto & o : obstacleBlocks) o.second->allocate_surface();
}

void IF3D_FishOperator::interpolateOnSkin(const double time, const int stepID, bool _dumpWake)
{
  #ifdef __useSkin_
  if(not(quaternion[1] == 0 && quaternion[2] == 0)){
    printf("the fish skin works only if the fish angular velocity is limited to the z axis. Aborting"); fflush(NULL);
    abort();
  }
  assert(quaternion[1] == 0 && quaternion[2] == 0);
  sr.updateStepId(stepID+obstacleID);
  myFish->computeSkinNormals(_2Dangle, CoM_interpolated);

  sr.nearestGridPoints(obstacleBlocks, vInfo, myFish->upperSkin->Npoints-1,
                myFish->upperSkin->midX,      myFish->upperSkin->midY,
                myFish->lowerSkin->midX,      myFish->lowerSkin->midY,
                myFish->upperSkin->normXSurf, myFish->upperSkin->normYSurf,
                myFish->lowerSkin->normXSurf, myFish->lowerSkin->normYSurf,
                position[2], vInfo[0].h_gridpoint, grid->getCartComm());

  #ifndef __RL_TRAINING
  //  if(rank==0) sr.print(obstacleID, stepID, time);
  #endif

  //if(_dumpWake && _uInf not_eq nullptr) dumpWake(stepID, time, _uInf);
  #endif
}

void IF3D_FishOperator::update(const int stepID, const double t, const double dt, const Real *Uinf)
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
  double P=2*(myFish->timeshift-myFish->time0/myFish->l_Tp)+myFish->phaseShift;
  sr.phaseShift = fmod(P,2)<0 ? 2+fmod(P,2) : fmod(P,2);
}

void IF3D_FishOperator::_parseArguments(ArgumentParser & parser)
{
  IF3D_ObstacleOperator::_parseArguments(parser);
  parser.set_strict_mode();
  parser.unset_strict_mode();
  Tperiod = parser("-T").asDouble(1.0);

  const bool optimizeT = parser("-optimizeT").asBool(false);
  if(optimizeT){

    char line[256];
    FILE *f = fopen("hingedParams.txt", "r");
    if(f==NULL){
      printf("hingedParams not found!\n"); abort();
    }

    double tvalue=0.0;
    int line_no = 0;

    while (fgets(line, 256, f)!= NULL) {
      if ((line[0] == '#')||(strlen(line)==0)) {
        printf("ignoring line %d\n", line_no);
        continue;
      }

      if (strstr(line, "T=")) {
        sscanf(line, "T=%lf", &tvalue);
      }
      line_no++;
    }
    if(tvalue==0.0){
      printf("Correct T=?? not found in hingedParams.txt\n");
      fflush(0);
      abort();
    }
    Tperiod = tvalue;
    printf("Tperiod overwritten: %f\n", Tperiod);
    assert(Tperiod >0.0);
    fclose(f);
  }

  phaseShift = parser("-phi").asDouble(0.0);

  //PID knobs
  bCorrectTrajectory = parser("-Correct").asBool(false);
  followX = parser("-followX").asDouble(-1);
  followY = parser("-followY").asDouble(-1);

  #ifdef __useSkin_
   bHasSkin = true;
  #endif
}

void IF3D_FishOperator::getCenterOfMass(double CM[3]) const
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

void IF3D_FishOperator::save(const int stepID, const double t, string filename)
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

void IF3D_FishOperator::restart(const double t, string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    if(!restartstream.good()){
      printf("Could not restart from file\n");
      return;
    }
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
