//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "Simulation.h"
#include <HDF5Dumper_MPI.h>
#include "ProcessOperatorsOMP.h"
#include <chrono>

void Simulation::_ic()
{
  CoordinatorIC coordIC(grid);
  profiler.push_start(coordIC.getName());
  coordIC(0);
  profiler.pop_stop();
}

void Simulation::setupGrid()
{
  parser.set_strict_mode();
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();
  parser.unset_strict_mode();
  nprocsy = parser("-nprocsy").asInt(1);
  parser.set_strict_mode();
  if(nprocsy==1) MPI_Comm_size(app_comm, &nprocsx);
  else           nprocsx = parser("-nprocsx").asInt();
  parser.unset_strict_mode();
  nprocsz = 1;

  if( not(bpdx%nprocsx==0 && bpdy%nprocsy==0 && bpdz%nprocsz==0) ){
  printf("Incompatible domain decomposition: bpd*/nproc* should be an integer");
  MPI_Abort(grid->getCartComm(), 1);
  }

  bpdx /= nprocsx; bpdy /= nprocsy; bpdz /= nprocsz;
  grid = new FluidGridMPI(nprocsx,nprocsy,nprocsz, bpdx,bpdy,bpdz, 1,app_comm);
  //dump = new  DumpGridMPI(nprocsx,nprocsy,nprocsz, bpdx,bpdy,bpdz, 1,app_comm);
  assert(grid != NULL);
  vInfo = grid->getBlocksInfo();

  MPI_Comm_rank(grid->getCartComm(),&rank);
  MPI_Comm_size(grid->getCartComm(),&nprocs);
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  const int nthreads = omp_get_max_threads();
  printf("Rank %d (of %d) with %d threads on host Hostname: %s\n",
          rank, nprocs, nthreads, hostname);
  //if (communicator not_eq nullptr) //Yo dawg I heard you like communicators.
  //  communicator->comm_MPI = grid->getCartComm();
  fflush(0);
  if(rank==0) {
    printf("Blocks per dimension: [%d %d %d]\n",bpdx,bpdy,bpdz);
    printf("Nranks per dimension: [%d %d %d]\n",nprocsx,nprocsy,nprocsz);
  }

  // setup 2D slices
  m_slices = SliceType::getSlices<SliceType>(parser, *grid);
}

void Simulation::parseArguments()
{
  nu = parser("-nu").asDouble();
  //length = parser("-length").asDouble();
  assert(nu>=0);
  parser.unset_strict_mode();
  bRestart = parser("-restart").asBool(false);
  b2Ddump = parser("-dump2D").asBool(false);
  b3Ddump = parser("-dump3D").asBool(true);
  DLM = parser("-use-dlm").asDouble(0);
  int dumpFreq = parser("-fdump").asDouble(0);  // dumpFreq==0 means dump freq (in #steps) is not active
  double dumpTime = parser("-tdump").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  saveTime = parser("-tsave").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if(saveFreq<=0 && dumpFreq>0) saveFreq = dumpFreq;
  if(saveTime<=0 && dumpTime>0) saveTime = dumpTime;

  nsteps = parser("-nsteps").asInt(0);    // nsteps==0   means stopping criteria is not active
  endTime = parser("-tend").asDouble(0);    // endTime==0  means stopping criteria is not active

  path2file = parser("-file").asString("./paternoster");
  path4serialization = parser("-serialization").asString("./");
  lambda = parser("-lambda").asDouble(1e6);

  CFL = parser("-CFL").asDouble(.1);
  uinf[0] = parser("-uinfx").asDouble(0.0);
  uinf[1] = parser("-uinfy").asDouble(0.0);
  uinf[2] = parser("-uinfz").asDouble(0.0);
  verbose = parser("-verbose").asBool(true);
}

void Simulation::setupObstacles()
{
  IF3D_ObstacleFactory obstacleFactory(grid, uinf);
  obstacle_vector = new IF3D_ObstacleVector(grid, obstacleFactory.create(parser));
  parser.unset_strict_mode();

  #ifdef RL_LAYER
    if(task not_eq nullptr) task->initializeObstacles(obstacle_vector);
  #endif

  if (rank == 0) {
    const double maxU = std::max({uinf[0], uinf[1], uinf[2]});
    const double length = obstacle_vector->getD();
    const double re = length * std::max(maxU, length) / nu;
    assert(length>0);
    printf("Kinematic viscosity: %f, Re: %f, length scale: %f\n",nu,re,length);
  }
}

void Simulation::setupOperators()
{
    pipeline.clear();
    pipeline.push_back(new CoordinatorComputeShape(grid, &obstacle_vector, &step, &time, uinf));
    pipeline.push_back(new CoordinatorPenalization(grid, &obstacle_vector, &lambda, uinf));
    pipeline.push_back(new CoordinatorComputeDiagnostics(grid, &obstacle_vector, &step, &time, &lambda, uinf));
    // For correct behavior Advection must always be followed Diffusion!
    pipeline.push_back(new CoordinatorAdvection<LabMPI>(uinf, grid));
    pipeline.push_back(new CoordinatorDiffusion<LabMPI>(nu, grid));
    pipeline.push_back(new CoordinatorPressure<LabMPI>(grid, &obstacle_vector));
    pipeline.push_back(new CoordinatorComputeForces(grid, &obstacle_vector, &step, &time, &nu, &bDump, uinf));
    //#ifndef _OPEN_BC_
    pipeline.push_back(new CoordinatorFadeOut(grid));
    //#endif
    if(rank==0) {
      cout << "Coordinator/Operator ordering:\n";
      for (int c=0; c<pipeline.size(); c++) cout << "\t" << pipeline[c]->getName() << endl;
    }
    //immediately call create!
    (*pipeline[0])(0);
}

double Simulation::calcMaxTimestep()
{
  double local_maxU = (double)findMaxUOMP(vInfo,*grid,uinf);
  double global_maxU;
  const double h = vInfo[0].h_gridpoint;

  MPI_Allreduce(&local_maxU, &global_maxU, 1, MPI::DOUBLE, MPI::MAX, grid->getCartComm());
  const double dtFourier = CFL*h*h/nu;
  const double dtCFL     = CFL*h/(std::fabs(global_maxU)+1e-8);
  double dt = std::min(dtCFL, dtFourier);

  // if DLM>=1, adapt lambda such that penal term is independent of time step
  // Probably best not used unless DLM>=10-100. Avoided during ramp-up (which
  // is the point of ramp-up: gradual insertion of obstacle)
  if (DLM >= 1) lambda = DLM / dt;

  if(!rank && verbose)
    printf("maxU %f dtF %f dtC %f dt %f\n", global_maxU,dtFourier,dtCFL,dt);

  if ( rampup && step<1000 ) {
    const double dt_max =  1e2*CFL*h;
    const double dt_min = 1e-2*CFL*h;
    //const double dt_ramp = dt_min + (dt_max-dt_min)*std::pow(step/1000.0, 2);
    const double dt_ramp = dt_min + (dt_max-dt_min)*(step/1000.0);
    if (dt_ramp<dt) {
      dt = dt_ramp;
      if(!rank && verbose)
        printf("Dt bounded by ramp-up: dt_ramp=%f\n",dt_ramp);
    }
  }
  return dt;
}

void Simulation::_serialize(const string append)
{
  if(!bDump) return;

  stringstream ssR;
  if (append == "") ssR<<"restart_";
  else ssR<<append;
  ssR<<std::setfill('0')<<std::setw(9)<<step;
  if (rank==0) cout<<"Saving to "<<path4serialization<<"/"<<ssR.str()<<endl;

  if (rank==0) { //rank 0 saves step id and obstacles
    obstacle_vector->save(step, time, path4serialization+"/"+ssR.str());
    //safety status in case of crash/timeout during grid save:
    string statusname = path4serialization+"/"+ssR.str()+".status";
    FILE * f = fopen(statusname.c_str(), "w");
    assert(f != NULL);
    fprintf(f, "time: %20.20e\n", time);
    fprintf(f, "stepid: %d\n", (int)step);
    fprintf(f, "uinfx: %20.20e\n", uinf[0]);
    fprintf(f, "uinfy: %20.20e\n", uinf[1]);
    fprintf(f, "uinfz: %20.20e\n", uinf[2]);
    fclose(f);
  }

  #if defined(_USE_HDF_)
    if(b2Ddump) {
      stringstream ssF;
      if (append == "")
       ssF<<"avemaria_"<<std::setfill('0')<<std::setw(9)<<step;
      else
       ssF<<"2D_"<<append<<std::setfill('0')<<std::setw(9)<<step;

      for (const auto& slice : m_slices)
        DumpSliceHDF5MPI<SliceType,StreamerVelocityVector>(slice, step, time, ssF.str(), path4serialization);
    }

    if(b3Ddump) {
      DumpHDF5_MPI<FluidGridMPI,StreamerVelocityVector>(*grid, step, time, ssR.str(), path4serialization);
      DumpHDF5_MPI<FluidGridMPI,StreamerChi>(*grid, step, time, ssR.str(), path4serialization);
    }
  #endif

  if (rank==0)
  { //saved the grid! Write status to remember most recent ping
    string restart_status = path4serialization+"/restart.status";
    FILE * f = fopen(restart_status.c_str(), "w");
    assert(f != NULL);
    fprintf(f, "time: %20.20e\n", time);
    fprintf(f, "stepid: %d\n", (int)step);
    fprintf(f, "uinfx: %20.20e\n", uinf[0]);
    fprintf(f, "uinfy: %20.20e\n", uinf[1]);
    fprintf(f, "uinfz: %20.20e\n", uinf[2]);
    fclose(f);
    printf("time:  %20.20e\n", time);
    printf("stepid: %d\n", (int)step);
    printf("uinfx: %20.20e\n", uinf[0]);
    printf("uinfy: %20.20e\n", uinf[1]);
    printf("uinfz: %20.20e\n", uinf[2]);
  }

  //CoordinatorDiagnostics coordDiags(grid,time,step);
  //coordDiags(dt);

  obstacle_vector->interpolateOnSkin(time, step);
}

void Simulation::_deserialize()
{
  {
    string restartfile = path4serialization+"/restart.status";
    FILE * f = fopen(restartfile.c_str(), "r");
    if (f == NULL) {
      printf("Could not restart... starting a new sim.\n");
      return;
    }
    assert(f != NULL);
    bool ret = true;
    ret = ret && 1==fscanf(f, "time: %e\n",   &time);
    ret = ret && 1==fscanf(f, "stepid: %d\n", &step);
    ret = ret && 1==fscanf(f, "uinfx: %e\n", &uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %e\n", &uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %e\n", &uinf[2]);
    fclose(f);
    if( (not ret) || step<0 || time<0) {
      printf("Error reading restart file. Aborting...\n");
      MPI_Abort(grid->getCartComm(), 1);
    }
  }

  stringstream ssR;
  ssR<<"restart_"<<std::setfill('0')<<std::setw(9)<<step;
  if (rank==0) cout << "Restarting from " << ssR.str() << endl;

  #if defined(_USE_HDF_)
    ReadHDF5_MPI<FluidGridMPI, StreamerVelocityVector>(*grid, ssR.str(),path4serialization);
  #else
    printf("Unable to restart without  HDF5 library. Aborting...\n");
    MPI_Abort(grid->getCartComm(), 1);
  #endif

  obstacle_vector->restart(time,path4serialization+"/"+ssR.str());

  printf("DESERIALIZATION: time is %f and step id is %d\n", time, (int)step);
  // prepare time for next save
  nextSaveTime = time + saveTime;
}

void Simulation::init()
{
  parseArguments();
  setupGrid();
  setupObstacles();

  if(bRestart) _deserialize();
  else _ic();

  setupOperators();

  MPI_Barrier(grid->getCartComm());
}

void Simulation::simulate()
{
    for (;;) {
        profiler.push_start("DT");
        const double dt = calcMaxTimestep();
        profiler.pop_stop();

        if (timestep(dt))
            break;
    }
}

bool Simulation::timestep(const double dt)
{
    bDump = (saveFreq>0 && (step+ 1)%saveFreq==0) ||
            (saveTime>0 && (time+dt)>nextSaveTime);
    if (bDump) nextSaveTime += saveTime;

    #ifdef RL_LAYER
      if(task not_eq nullptr)
        if((*task)(step, time)) {
          if(rank==0)
          cout<<"Finished RL task at time "<<time<<endl;
          return;
        }
    #endif

    for (int c=0; c<pipeline.size(); c++) {
      profiler.push_start(pipeline[c]->getName());
      (*pipeline[c])(dt);
      profiler.pop_stop();
    }
    step++;
    time+=dt;

    if(!rank && verbose)
      printf("%d : %f uInf {%f %f %f}\n",step,time,uinf[0],uinf[1],uinf[2]);

    profiler.push_start("Save");
    _serialize();
    profiler.pop_stop();

    if (step % 50 == 0 && !rank && verbose) profiler.printSummary();
    if ((endTime>0 && time>endTime) || (nsteps!=0 && step>=nsteps))
    {
      if(rank==0)
      cout<<"Finished at time "<<time<<" in "<<step<<" step of "<<nsteps<<endl;
      return true;  // Finished.
    }

    return false;  // Not yet finished.
}
