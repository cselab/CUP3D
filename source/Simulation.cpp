//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//
#include "Simulation.h"

#include "operators/InitialConditions.h"
#include "operators/ObstaclesCreate.h"
#include "operators/AdvectionDiffusion.h"
#include "operators/ObstaclesUpdate.h"
#include "operators/Penalization.h"
#include "operators/PressureRHS.h"
#include "operators/PressureProjection.h"
#include "operators/ComputeDissipation.h"
#include "operators/FluidSolidForces.h"
#include "operators/ProcessHelpers.h"

#include "obstacles/ObstacleVector.h"
#include "obstacles/ObstacleFactory.h"

#include <Cubism/HDF5Dumper_MPI.h>
#include <Cubism/ArgumentParser.h>

#include <iomanip>
#include <iostream>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

Simulation::Simulation(MPI_Comm mpicomm, ArgumentParser & parser) : sim(mpicomm, parser) {
  // Make sure given arguments are valid
  if( sim.verbose )
    std::cout << "[CUP3D] Parsing Arguments.. " << std::endl;
  sim._preprocessArguments();

  // Setup and Initialize Grid
  if( sim.verbose )
    std::cout << "[CUP3D] Allocating Grid.. " << std::endl;
  setupGrid();
  touch();

  // Setup Computational Pipeline
  if( sim.verbose )
    std::cout << "[CUP3D] Creating Computational Pipeline.. " << std::endl;
  setupOperators();

  // Initalize Obstacles
  if( sim.verbose )
    std::cout << "[CUP3D] Initializing Obstacles.. " << std::endl;
  sim.obstacle_vector = new ObstacleVector(sim);
  ObstacleFactory(sim).addObstacles(parser);

  // CreateObstacles
  if( sim.verbose )
    std::cout << "[CUP3D] Creating Obstacles.. " << std::endl;
  (*sim.pipeline[0])(0);

  // Save Initial Flow Field to File
  _serialize("init");

  // Initialize Flow Field
  if( sim.verbose )
    std::cout << "[CUP3D] Initializing Flow Field.. " << std::endl;
  const bool bRestart = parser("-restart").asBool(false);
  if (bRestart)
    _deserialize();
  else if (sim.icFromH5 != "")
    _icFromH5(sim.icFromH5);
  else
    _ic();

  // Initial refinement of Grid
  if( sim.verbose )
    std::cout << "[CUP3D] Performing Initial Refinement of Grid.. " << std::endl;
  refineGrid();
}

void Simulation::refineGrid()
{
  // Initial Compression of Grid
  for (int l = 0 ; l < sim.levelMax ; l++)
  {
    if( sim.verbose )
      std::cout << "[CUP3D] - refinement " << l << "/" << sim.levelMax-1 << std::endl;
    // CreateObstacles
    (*sim.pipeline[0])(0);

    // Compression of Grid
    sim.amr->Tag();
    sim.amr2->TagLike(sim.vInfo());
    sim.amr->Adapt(sim.time,sim.verbose,false);
    sim.amr2->Adapt(sim.time,false,true);

    // After mesh is refined/compressed min_pos and max_pos need to change
    const std::vector<BlockInfo>& vInfo = sim.vInfo();
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      b.min_pos = vInfo[i].pos<Real>(0, 0, 0);
      b.max_pos = vInfo[i].pos<Real>(FluidBlock::sizeX-1,FluidBlock::sizeY-1,FluidBlock::sizeZ-1);
    }

    //This may not be needed but has zero cost 
    if (l != sim.levelMax-1) touch();

  }

  // Save Initial Flow Field to File
  _serialize("refined");
}

const std::vector<std::shared_ptr<Obstacle>>& Simulation::getObstacleVector() const
{
    return sim.obstacle_vector->getObstacleVector();
}

void Simulation::_ic()
{
  InitialConditions coordIC(sim);
  sim.startProfiler(coordIC.getName());
  coordIC(0);
  sim.stopProfiler();
}

void Simulation::_icFromH5(std::string h5File)
{
  if (sim.rank==0) std::cout << "Extracting Initial Conditions from " << h5File << std::endl;

  #ifdef CUBISM_USE_HDF
    ReadHDF5_MPI<StreamerVelocityVector, DumpReal>(* sim.grid,
      h5File, sim.path4serialization);
  #else
    printf("Unable to restart without  HDF5 library. Aborting...\n");
    fflush(0); MPI_Abort(sim.grid->getCartComm(), 1);
  #endif

  sim.obstacle_vector->restart(sim.path4serialization+"/"+sim.icFromH5);

  // prepare time for next save
  sim.nextSaveTime = sim.time + sim.saveTime;
  MPI_Barrier(sim.app_comm);
}

void Simulation::setupGrid()
{
  // if(sim.rank==0) printf("Grid of sizes: %f %f %f\n", sim.extent[0],sim.extent[1],sim.extent[2]);
  sim.grid = new FluidGridMPI(1, //these arguments are not used in Cubism-AMR
                              1, //these arguments are not used in Cubism-AMR
                              1, //these arguments are not used in Cubism-AMR
                              sim.bpdx,
                              sim.bpdy,
                              sim.bpdz,
                              sim.maxextent,
                              sim.levelStart,sim.levelMax,sim.app_comm,
                              (sim.BCx_flag == periodic),
                              (sim.BCy_flag == periodic),
                              (sim.BCz_flag == periodic));
  assert(sim.grid != nullptr);

  sim.gridPoisson = new FluidGridMPIPoisson(1, //these arguments are not used in Cubism-AMR
                                            1, //these arguments are not used in Cubism-AMR
                                            1, //these arguments are not used in Cubism-AMR
                                            sim.bpdx,
                                            sim.bpdy,
                                            sim.bpdz,
                                            sim.maxextent,
                                            sim.levelStart,sim.levelMax,sim.app_comm,
                                            (sim.BCx_flag == periodic),
                                            (sim.BCy_flag == periodic),
                                            (sim.BCz_flag == periodic));

  //Refine/compress only according to chi field for now
  sim.amr = new AMR( *(sim.grid),sim.Rtol,sim.Ctol);
  sim.amr2 = new AMR2( *(sim.gridPoisson),sim.Rtol,sim.Ctol);

  const std::vector<BlockInfo>& vInfo = sim.vInfo();
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<vInfo.size(); i++) {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    b.min_pos = vInfo[i].pos<Real>(0, 0, 0);
    b.max_pos = vInfo[i].pos<Real>(FluidBlock::sizeX-1,
                                   FluidBlock::sizeY-1,
                                   FluidBlock::sizeZ-1);
  }
}

void Simulation::setupOperators()
{
  // Creates the char function, sdf, and def vel for all obstacles at the curr
  // timestep. At this point we do NOT know the translation and rot vel of the
  // obstacles. We need to solve implicit system when the pre-penalization vel
  // is finalized on the grid.
  // Here we also compute obstacles' centres of mass which are computed from
  // the char func on the grid. This is different from "position" which is
  // the quantity that is advected and used to construct shape.
  sim.pipeline.push_back(new CreateObstacles(sim));

  // Performs:
  // \tilde{u} = u_t + \delta t (\nu \nabla^2 u_t - (u_t \cdot \nabla) u_t )
  sim.pipeline.push_back(new AdvectionDiffusion(sim));

  // Update obstacle velocities and penalize velocity
  sim.pipeline.push_back(new UpdateObstacles(sim));
  sim.pipeline.push_back(new Penalization(sim));

  // Places Udef on the grid and computes the RHS of the Poisson Eq
  // overwrites tmpU, tmpV, tmpW and pressure solver's RHS
  sim.pipeline.push_back(new PressureRHS(sim));

  // Solves the Poisson Eq to get the pressure and finalizes the velocity
  sim.pipeline.push_back(new PressureProjection(sim));

  // With finalized velocity and pressure, compute forces and dissipation
  sim.pipeline.push_back(new ComputeForces(sim));
  sim.pipeline.push_back(new ComputeDissipation(sim));

  //sim.pipeline.push_back(new ComputeDivergence(sim));
  if(sim.rank==0) {
    printf("[CUP3D] Operator ordering:\n");
    for (size_t c=0; c<sim.pipeline.size(); c++)
      printf("\t - %s\n", sim.pipeline[c]->getName().c_str());
  }
}

void Simulation::_serialize(const std::string append)
{
  sim.startProfiler("DumpHDF5_MPI");

  sim.grid->UpdateMyGroups();

  std::stringstream name;
  if (append == "") name<<"restart_";
  else name<<append;
  name<<std::setfill('0')<<std::setw(9)<<sim.step;
  auto * grid2Dump = sim.grid;
  if (sim.dumpOmega || sim.dumpOmegaX || sim.dumpOmegaY || sim.dumpOmegaZ)
  {
    ComputeVorticity  FindOmega(sim);
    FindOmega(0);
  }
  if (sim.dumpP        ) DumpHDF5_MPI<StreamerPressure      , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerPressure       ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpChi      ) DumpHDF5_MPI<StreamerChi           , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerChi            ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmega    ) DumpHDF5_MPI<StreamerTmpVector     , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVector      ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaX   ) DumpHDF5_MPI<StreamerTmpVectorX    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorX     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaY   ) DumpHDF5_MPI<StreamerTmpVectorY    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorY     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaZ   ) DumpHDF5_MPI<StreamerTmpVectorZ    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorZ     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocity ) DumpHDF5_MPI<StreamerVelocityVector, DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelocityVector ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityX) DumpHDF5_MPI<StreamerVelVectorX    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorX     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityY) DumpHDF5_MPI<StreamerVelVectorY    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorY     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityZ) DumpHDF5_MPI<StreamerVelVectorZ    , DumpReal, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorZ     ::prefix() + name.str(),sim.path4serialization);

  sim.stopProfiler();
}

void Simulation::_deserialize()
{
  {
    std::string restartfile = sim.path4serialization+"/restart.status";
    FILE * f = fopen(restartfile.c_str(), "r");
    if (f == NULL) {
      printf("Could not restart... starting a new sim.\n");
      return;
    }
    assert(f != NULL);
    bool ret = true;
    ret = ret && 1==fscanf(f, "time: %le\n",   &sim.time);
    ret = ret && 1==fscanf(f, "stepid: %d\n", &sim.step);
    #ifndef CUP_SINGLE_PRECISION
    ret = ret && 1==fscanf(f, "uinfx: %le\n", &sim.uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %le\n", &sim.uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %le\n", &sim.uinf[2]);
    #else // CUP_SINGLE_PRECISION
    ret = ret && 1==fscanf(f, "uinfx: %e\n", &sim.uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %e\n", &sim.uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %e\n", &sim.uinf[2]);
    #endif // CUP_SINGLE_PRECISION
    fclose(f);
    if( (not ret) || sim.step<0 || sim.time<0) {
      printf("Error reading restart file. Aborting...\n");
      fflush(0); MPI_Abort(sim.grid->getCartComm(), 1);
    }
  }

  std::stringstream ssR;
  ssR<<"restart_"<<std::setfill('0')<<std::setw(9)<<sim.step;
  if (sim.rank==0) std::cout << "Restarting from " << ssR.str() << "\n";

  #ifdef CUBISM_USE_HDF
    ReadHDF5_MPI<StreamerVelocityVector, DumpReal>(* sim.grid,
      StreamerVelocityVector::prefix()+ssR.str(), sim.path4serialization);
  #else
    printf("Unable to restart without  HDF5 library. Aborting...\n");
    fflush(0); MPI_Abort(sim.grid->getCartComm(), 1);
  #endif

  sim.obstacle_vector->restart(sim.path4serialization+"/"+ssR.str());

  printf("DESERIALIZATION: time is %f and step id is %d\n", sim.time, sim.step);
  // prepare time for next save
  sim.nextSaveTime = sim.time + sim.saveTime;
}

void Simulation::run()
{
  for (;;) {
    sim.startProfiler("DT");
    const double dt = calcMaxTimestep();
    sim.stopProfiler();

    if (timestep(dt)) break;
  }
}

double Simulation::calcMaxTimestep()
{
  const double dt_old = sim.dt;
  sim.dt_old = sim.dt;
  const double hMin = sim.hmin;
  double CFL = sim.CFL;
  sim.uMax_measured = findMaxU(sim);

  if( CFL > 0 )
  {
    const double dtDiffusion = (1./ 6.) * ( hMin * hMin / sim.nu );
    const double dtAdvection = hMin / ( sim.uMax_measured + 1e-8 );
    if ( sim.step < sim.rampup )
    {
      const double x = sim.step / (double) sim.rampup;
      const double rampCFL = std::exp(std::log(1e-3)*(1-x) + std::log(CFL)*x);
      sim.dt = std::min(dtDiffusion, rampCFL * dtAdvection);
    }
    else
      sim.dt = std::min(dtDiffusion, CFL * dtAdvection);
  }
  else
  {
    CFL = ( sim.uMax_measured + 1e-8 ) * sim.dt / hMin;
  }

  if( sim.dt <= 0 ){
    fprintf(stderr, "dt <= 0. CFL=%f, hMin=%f, sim.uMax_measured=%f. Aborting...\n", CFL, hMin, sim.uMax_measured);
    fflush(0); MPI_Abort(sim.grid->getCartComm(), 1);
  }


  // if DLM>0, adapt lambda such that penal term is independent of time step
  if (sim.DLM > 0) sim.lambda = sim.DLM / sim.dt;

  if( sim.rank == 0 ) {
    printf("[CUP3D] step: %d, time: %f, dt: %.2e, uinf: {%f %f %f}, maxU:%f, minH:%f, CFL:%.2e, lambda:%.2e\n",sim.step,sim.time, sim.dt, sim.uinf[0],sim.uinf[1],sim.uinf[2], sim.uMax_measured, hMin, CFL, sim.lambda);
  }

  if (sim.step > sim.step_2nd_start)
  {
    const double a = dt_old;
    const double b = sim.dt;
    const double c1 = -(a+b)/(a*b);
    const double c2 = b/(a+b)/a;
    sim.coefU[0] = -b*(c1+c2);
    sim.coefU[1] = b*c1;
    sim.coefU[2] = b*c2;
    //sim.coefU[0] = 1.5;
    //sim.coefU[1] = -2.0;
    //sim.coefU[2] = 0.5;
  }
  return sim.dt;
}

bool Simulation::timestep(const double dt)
{
  const bool bDumpFreq = (sim.saveFreq>0 && (sim.step+ 1)%sim.saveFreq==0);
  const bool bDumpTime = (sim.saveTime>0 && (sim.time+dt)>sim.nextSaveTime);
  if (bDumpTime) sim.nextSaveTime += sim.saveTime;
  sim.bDump = (bDumpFreq || bDumpTime);

  const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  for (size_t c=0; c< sim.pipeline.size() ; c++)
  {
    (*sim.pipeline[c])(dt);
  }
  if (sim.step % 5 == 0 || sim.step < 10)
  {
    sim.amr->Tag();
    sim.amr2->TagLike(sim.vInfo());
    sim.amr->Adapt(sim.time,sim.verbose,false);
    sim.amr2->Adapt(sim.time,false,true);
    //After mesh is refined/coarsened the arrays min_pos and max_pos need to change.
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      b.min_pos = vInfo[i].pos<Real>(0, 0, 0);
      b.max_pos = vInfo[i].pos<Real>(FluidBlock::sizeX-1,FluidBlock::sizeY-1,FluidBlock::sizeZ-1);
    }
  }
  sim.step++;
  sim.time += dt;

  sim.startProfiler("Save");
  if( sim.bDump ) _serialize();
  sim.stopProfiler();

  if (sim.step % 50 == 0 && sim.verbose) sim.printResetProfiler();
  if ((sim.endTime>0 && sim.time>sim.endTime) ||
      (sim.nsteps!=0 && sim.step>=sim.nsteps) ) {
    if(sim.verbose)
    {
      sim.printResetProfiler();
      std::cout<<"Finished at time "<<sim.time<<" in "<<sim.step<<" steps.\n";
    }
    return true;  // Finished.
  }

  return false;  // Not yet finished.
}

CubismUP_3D_NAMESPACE_END
