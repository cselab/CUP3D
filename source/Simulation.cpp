//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//
#include "Simulation.h"

#include "Cubism/ArgumentParser.h"
#include "Cubism/HDF5Dumper_MPI.h"
#include "CoordinatorAdvectDiffuse.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorComputeDissipation.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorFadeOut.h"
#include "CoordinatorIC.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorPressure.h"
//#include "CoordinatorVorticity.h"
#include "IF3D_ObstacleFactory.h"
#include "ProcessOperatorsOMP.h"

/*
 * Initialization from cmdline arguments is done in few steps, because grid has
 * to be created before the obstacles and slices are created.
 */
Simulation::Simulation(MPI_Comm mpicomm, ArgumentParser &parser)
    : app_comm(mpicomm)
{
  MPI_Comm_rank(app_comm, &rank);
  if (rank == 0)
    parser.print_args();

  // ========== SIMULATION ==========
  // GRID
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();
  nprocsx = parser("-nprocsx").asInt(-1);
  nprocsy = parser("-nprocsy").asInt(-1);

  // FLOW
  nu = parser("-nu").asDouble();
  lambda = parser("-lambda").asDouble(1e6);
  CFL = parser("-CFL").asDouble(.1);
  uinf[0] = parser("-uinfx").asDouble(0.0);
  uinf[1] = parser("-uinfy").asDouble(0.0);
  uinf[2] = parser("-uinfz").asDouble(0.0);

  // PIPELINE
  computeDissipation = (bool)parser("-compute-dissipation").asInt(0);
#ifndef _UNBOUNDED_FFT_
  fadeOutLength = parser("-fade_len").asDouble(.005);
#endif

  // OUTPUT
  verbose = parser("-verbose").asBool(true);
  b2Ddump = parser("-dump2D").asBool(false);
  b3Ddump = parser("-dump3D").asBool(true);
  DLM = parser("-use-dlm").asDouble(0);

  int dumpFreq = parser("-fdump").asDouble(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  double dumpTime = parser("-tdump").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);         // dumpFreq==0 means dump freq (in #steps) is not active
  saveTime = parser("-tsave").asDouble(0.0);    // dumpTime==0 means dump freq (in time)   is not active

  nsteps = parser("-nsteps").asInt(0);    // 0 to disable this stopping critera.
  endTime = parser("-tend").asDouble(0);  // 0 to disable this stopping critera.

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if (saveFreq <= 0 && dumpFreq > 0) saveFreq = dumpFreq;
  if (saveTime <= 0 && dumpTime > 0) saveTime = dumpTime;

  path4serialization = parser("-serialization").asString("./");


  // ========= SETUP GRID ==========
  // Grid has to be initialized before slices and obstacles.
  _argumentsSanityCheck();
  setupGrid();

  // =========== SLICES ============
  #ifdef DUMPGRID
    m_slices = SliceType::getEntities<SliceType>(parser, *dump);
  #else
    m_slices = SliceType::getEntities<SliceType>(parser, *grid);
  #endif

  // ========== OBSTACLES ==========
  IF3D_ObstacleFactory factory(grid, uinf);
  setObstacleVector(new IF3D_ObstacleVector(grid, factory.create(parser)));

  // ====== SETUP OPERATORS ========
  setupOperators();

  bool restart = parser("-restart").asBool(false);
  if (restart)
    _deserialize();
  else
    _ic();

  MPI_Barrier(app_comm);
}

void Simulation::_ic()
{
  CoordinatorIC coordIC(grid);
  profiler.push_start(coordIC.getName());
  coordIC(0);
  profiler.pop_stop();
}

void Simulation::setupGrid()
{
  assert(bpdx > 0);
  assert(bpdy > 0);
  assert(bpdz > 0);
  MPI_Comm_rank(app_comm, &rank);
  MPI_Comm_size(app_comm, &nprocs);

  if (nprocsy < 0)
    nprocsy = 1;
  if (nprocsy == 1) {
    nprocsx = nprocs;  // Override existing value!
  } else if (nprocsx < 0) {
    nprocsx = nprocs / nprocsy;
  }
  nprocsz = 1;

  if (nprocsx * nprocsy * nprocsz != nprocs) {
    fprintf(stderr, "Invalid domain decomposition. %d x %d x %d != %d!\n",
            nprocsx, nprocsy, nprocsz, nprocs);
    MPI_Abort(app_comm, 1);
  }

  if (bpdx % nprocsx != 0 || bpdy % nprocsy != 0 || bpdz % nprocsz !=0) {
    printf("Incompatible domain decomposition: bpd*/nproc* should be an integer");
    MPI_Abort(app_comm, 1);
  }

  bpdx /= nprocsx;
  bpdy /= nprocsy;
  bpdz /= nprocsz;
  grid = new FluidGridMPI(nprocsx,nprocsy,nprocsz, bpdx,bpdy,bpdz, 1,app_comm);
  assert(grid != NULL);
  vInfo = grid->getBlocksInfo();

  #ifdef DUMPGRID
    // create new comm so that if there is a barrier main work is not affected
    MPI_Comm_split(app_comm, 0, rank, &dump_comm);
    dump = new  DumpGridMPI(nprocsx,nprocsy,nprocsz, bpdx,bpdy,bpdz, 1, dump_comm);
  #endif

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

}

void Simulation::_argumentsSanityCheck()
{
  // Flow.
  assert(nu >= 0);
  assert(lambda > 0.0);
  assert(CFL > 0.0);

  // Output.
  assert(saveFreq >= 0.0);
  assert(saveTime >= 0.0);
}

void Simulation::setObstacleVector(IF3D_ObstacleVector * const obstacle_vector_)
{
  obstacle_vector = obstacle_vector_;
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

  // For correct behavior Advection must always precede Diffusion!
  // pipeline.push_back(new CoordinatorAdvection<LabMPI>(uinf, grid));
  // pipeline.push_back(new CoordinatorDiffusion<LabMPI>(nu, grid));
  pipeline.push_back(new CoordinatorAdvectDiffuse<LabMPI>(nu, uinf, grid));

  pipeline.push_back(new CoordinatorPressure<LabMPI>(grid, &obstacle_vector));
  pipeline.push_back(new CoordinatorComputeForces(grid, &obstacle_vector, &step, &time, &nu, &bDump, uinf));
  if(computeDissipation)
    pipeline.push_back(new CoordinatorComputeDissipation<LabMPI>(grid,nu,&step,&time));

  #ifndef _UNBOUNDED_FFT_
    pipeline.push_back(new CoordinatorFadeOut(grid, fadeOutLength));
  #endif /* _UNBOUNDED_FFT_ */

  if(rank==0) {
    printf("Coordinator/Operator ordering:\n");
    for (size_t c=0; c<pipeline.size(); c++)
      printf("\t%s\n", pipeline[c]->getName().c_str());
  }
  //immediately call create!
  (*pipeline[0])(0);
}

double Simulation::calcMaxTimestep()
{
  double local_maxU = (double)findMaxUOMP(vInfo,*grid,uinf);
  double global_maxU;
  const double h = vInfo[0].h_gridpoint;

  MPI_Allreduce(&local_maxU, &global_maxU, 1, MPI_DOUBLE, MPI_MAX, grid->getCartComm());
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

void Simulation::_serialize(const std::string append)
{
  if(!bDump) return;

  std::stringstream ssR;
  if (append == "") ssR<<"restart_";
  else ssR<<append;
  ssR<<std::setfill('0')<<std::setw(9)<<step;
  if(rank==0) std::cout<<"Saving to "<<path4serialization<<"/"<<ssR.str()<<"\n";

  if (rank==0) { //rank 0 saves step id and obstacles
    obstacle_vector->save(step, time, path4serialization+"/"+ssR.str());
    //safety status in case of crash/timeout during grid save:
    std::string statusname = path4serialization+"/"+ssR.str()+".status";
    FILE * f = fopen(statusname.c_str(), "w");
    assert(f != NULL);
    fprintf(f, "time: %20.20e\n", time);
    fprintf(f, "stepid: %d\n", (int)step);
    fprintf(f, "uinfx: %20.20e\n", uinf[0]);
    fprintf(f, "uinfy: %20.20e\n", uinf[1]);
    fprintf(f, "uinfz: %20.20e\n", uinf[2]);
    fclose(f);
  }

  #ifdef CUBISM_USE_HDF
  std::stringstream ssF;
  if (append == "")
   ssF<<"avemaria_"<<std::setfill('0')<<std::setw(9)<<step;
  else
  ssF<<"2D_"<<append<<std::setfill('0')<<std::setw(9)<<step;

  #ifdef DUMPGRID
    // if a thread was already created, make sure it has finished
    if(dumper not_eq nullptr) {
      dumper->join();
      delete dumper;
      dumper = nullptr;
    }
    // copy qois from grid to dump
    copyDumpGrid(*grid, *dump);
    const auto & grid2Dump = * dump;
  #else //DUMPGRID
    const auto & grid2Dump = * grid;
  #endif //DUMPGRID

  const auto name3d = ssR.str(), name2d = ssF.str(); // sstreams are weird

  const auto dumpFunction = [=] () {
    if(b2Ddump) {
      for (const auto& slice : m_slices) {
        DumpSliceHDF5MPI<StreamerVelocityVector, DumpReal>(
            slice, step, time, StreamerVelocityVector::prefix()+name2d, path4serialization);
        DumpSliceHDF5MPI<StreamerPressure, DumpReal>(
            slice, step, time, StreamerPressure::prefix()+name2d, path4serialization);
        DumpSliceHDF5MPI<StreamerChi, DumpReal>(
            slice, step, time, StreamerChi::prefix()+name2d, path4serialization);
      }
    }
    if(b3Ddump) {
      DumpHDF5_MPI<StreamerVelocityVector, DumpReal>(
          grid2Dump, step, time, StreamerVelocityVector::prefix()+name3d, path4serialization);
      DumpHDF5_MPI<StreamerPressure, DumpReal>(
          grid2Dump, step, time, StreamerPressure::prefix()+name3d, path4serialization);
      DumpHDF5_MPI<StreamerChi, DumpReal>(
          grid2Dump, step, time, StreamerChi::prefix()+name3d, path4serialization);
    }
  };

  #ifdef DUMPGRID
    dumper = new std::thread( dumpFunction );
  #else //DUMPGRID
    dumpFunction();
  #endif //DUMPGRID
  #endif //CUBISM_USE_HDF


  if (rank==0)
  { //saved the grid! Write status to remember most recent ping
    std::string restart_status = path4serialization+"/restart.status";
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
  #ifdef RL_LAYER
  obstacle_vector->interpolateOnSkin(time, step);
  #endif
}

void Simulation::_deserialize()
{
  {
    std::string restartfile = path4serialization+"/restart.status";
    FILE * f = fopen(restartfile.c_str(), "r");
    if (f == NULL) {
      printf("Could not restart... starting a new sim.\n");
      return;
    }
    assert(f != NULL);
    bool ret = true;
    ret = ret && 1==fscanf(f, "time: %le\n",   &time);
    ret = ret && 1==fscanf(f, "stepid: %d\n", &step);
    #ifndef _FLOAT_PRECISION_
    ret = ret && 1==fscanf(f, "uinfx: %le\n", &uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %le\n", &uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %le\n", &uinf[2]);
    #else // _FLOAT_PRECISION_
    ret = ret && 1==fscanf(f, "uinfx: %e\n", &uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %e\n", &uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %e\n", &uinf[2]);
    #endif // _FLOAT_PRECISION_
    fclose(f);
    if( (not ret) || step<0 || time<0) {
      printf("Error reading restart file. Aborting...\n");
      MPI_Abort(grid->getCartComm(), 1);
    }
  }

  std::stringstream ssR;
  ssR<<"restart_"<<std::setfill('0')<<std::setw(9)<<step;
  if (rank==0) std::cout << "Restarting from " << ssR.str() << "\n";

  #ifdef CUBISM_USE_HDF
    ReadHDF5_MPI<StreamerVelocityVector, DumpReal>(*grid,
      StreamerVelocityVector::prefix()+ssR.str(), path4serialization);
  #else
    printf("Unable to restart without  HDF5 library. Aborting...\n");
    MPI_Abort(grid->getCartComm(), 1);
  #endif

  obstacle_vector->restart(time, path4serialization+"/"+ssR.str());

  printf("DESERIALIZATION: time is %f and step id is %d\n", time, (int)step);
  // prepare time for next save
  nextSaveTime = time + saveTime;
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
    const bool bDumpFreq = (saveFreq>0 && (step+ 1)%saveFreq==0);
    const bool bDumpTime = (saveTime>0 && (time+dt)>nextSaveTime);
    if (bDumpTime) nextSaveTime += saveTime;
    bDump = (bDumpFreq || bDumpTime);

    #ifdef RL_LAYER
      if(task not_eq nullptr)
        if((*task)(step, time)) {
          if(rank==0)
          cout<<"Finished RL task at time "<<time<<endl;
          return;
        }
    #endif

    for (size_t c=0; c<pipeline.size(); c++) {
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

    if (step % 10 == 0 && !rank && verbose) {
      profiler.printSummary();
      profiler.reset();
    }
    if ((endTime>0 && time>endTime) || (nsteps!=0 && step>=nsteps))
    {
      if(rank==0)
      std::cout<<"Finished at time "<<time<<" in "<<step<<" steps.\n";
      return true;  // Finished.
    }

    return false;  // Not yet finished.
}
