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

#include "Obstacles/ObstacleVector.h"
#include "Obstacles/ObstacleFactory.h"

#include <Cubism/HDF5Dumper_MPI.h>
#include <Cubism/ArgumentParser.h>

#include <iomanip>
#include <iostream>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

std::shared_ptr<Simulation> createSimulation(
    const MPI_Comm comm,
    const std::vector<std::string> &argv)
{
  std::vector<char *> cargv(argv.size() + 1);
  char cmd[] = "prg";
  cargv[0] = cmd;
  for (size_t i = 0; i < argv.size(); ++i)
    cargv[i + 1] = const_cast<char *>(argv[i].data());
  ArgumentParser parser((int)cargv.size(), cargv.data());
  return std::make_shared<Simulation>(comm, parser);
}

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
  initialGridRefinement();
}

void Simulation::initialGridRefinement()
{
  for (int l = 0 ; l < 3*sim.levelMax ; l++)
  {
    if( sim.verbose )
      std::cout << "[CUP3D] - refinement " << l << "/" << 3*sim.levelMax-1 << std::endl;
    // CreateObstacles
    (*sim.pipeline[0])(0);

    // Refinement or compression of Grid
    adaptMesh();

    //This may not be needed but has zero cost
    if (l != 3*sim.levelMax-1) touch();
  }

  // Save Initial Flow Field to File
  //if ( sim.saveFreq>0 || sim.dumpTime>0 ) _serialize("init");
}

void Simulation::adaptMesh()
{
  sim.startProfiler("Mesh refinement");
  sim.amr->Tag();
  sim.lhs_amr->TagLike(sim.vInfo());
  sim.z_amr  ->TagLike(sim.vInfo());
  sim.amr->Adapt(sim.time,sim.verbose,false);
  sim.lhs_amr->Adapt(sim.time,false,true);
  sim.z_amr  ->Adapt(sim.time,false,true);
  sim.stopProfiler();
}

const std::vector<std::shared_ptr<Obstacle>>& Simulation::getShapes() const
{
    return sim.obstacle_vector->getObstacleVector();
}

void Simulation::_ic()
{
  InitialConditions coordIC(sim);
  coordIC(0);
}

void Simulation::_icFromH5(std::string h5File)
{
  if (sim.rank==0) std::cout << "Extracting Initial Conditions from " << h5File << std::endl;

  ReadHDF5_MPI<StreamerVelocityVector, Real>(* sim.grid, h5File, sim.path4serialization);

  sim.obstacle_vector->restart(sim.path4serialization+"/"+sim.icFromH5);

  // prepare time for next save
  sim.nextSaveTime = sim.time + sim.dumpTime;
  MPI_Barrier(sim.app_comm);
}

void Simulation::setupGrid()
{
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

  sim.lhs         = new ScalarGrid         (1, //these arguments are not used in Cubism-AMR
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

  sim.z           = new ScalarGrid         (1, //these arguments are not used in Cubism-AMR
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
  sim.lhs_amr = new ScalarAMR( *(sim.lhs),sim.Rtol,sim.Ctol);
  sim.z_amr   = new ScalarAMR( *(sim.z  ),sim.Rtol,sim.Ctol);
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
  sim.pipeline.push_back(std::make_shared<CreateObstacles>(sim));

  // Performs:
  // \tilde{u} = u_t + \delta t (\nu \nabla^2 u_t - (u_t \cdot \nabla) u_t )
  sim.pipeline.push_back(std::make_shared<AdvectionDiffusion>(sim));

  // Update obstacle velocities and penalize velocity
  sim.pipeline.push_back(std::make_shared<UpdateObstacles>(sim));
  sim.pipeline.push_back(std::make_shared<Penalization>(sim));

  // Places Udef on the grid and computes the RHS of the Poisson Eq
  // overwrites tmpU, tmpV, tmpW and pressure solver's RHS
  sim.pipeline.push_back(std::make_shared<PressureRHS>(sim));

  // Solves the Poisson Eq to get the pressure and finalizes the velocity
  sim.pipeline.push_back(std::make_shared<PressureProjection>(sim));

  // With finalized velocity and pressure, compute forces and dissipation
  sim.pipeline.push_back(std::make_shared<ComputeForces>(sim));
  sim.pipeline.push_back(std::make_shared<ComputeDissipation>(sim));

  //sim.pipeline.push_back(std::make_shared<ComputeDivergence>(sim));
  if(sim.rank==0) {
    printf("[CUP3D] Operator ordering:\n");
    for (size_t c=0; c<sim.pipeline.size(); c++)
      printf("\t - %s\n", sim.pipeline[c]->getName().c_str());
  }
}

void Simulation::_serialize(const std::string append)
{
  sim.startProfiler("DumpHDF5_MPI");

  std::stringstream name;
  if (append == "") name<<"restart_";
  else name<<append;
  name<<std::setfill('0')<<std::setw(9)<<sim.step;
  auto * grid2Dump = sim.grid;
  if (sim.dumpOmega || sim.dumpOmegaX || sim.dumpOmegaY || sim.dumpOmegaZ)
    computeVorticity();
  if (sim.dumpP        ) DumpHDF5_MPI2<StreamerPressure      , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerPressure       ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpChi      ) DumpHDF5_MPI2<StreamerChi           , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerChi            ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmega    ) DumpHDF5_MPI2<StreamerTmpVector     , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVector      ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaX   ) DumpHDF5_MPI2<StreamerTmpVectorX    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorX     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaY   ) DumpHDF5_MPI2<StreamerTmpVectorY    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorY     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpOmegaZ   ) DumpHDF5_MPI2<StreamerTmpVectorZ    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerTmpVectorZ     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocity ) DumpHDF5_MPI2<StreamerVelocityVector, Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelocityVector ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityX) DumpHDF5_MPI2<StreamerVelVectorX    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorX     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityY) DumpHDF5_MPI2<StreamerVelVectorY    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorY     ::prefix() + name.str(),sim.path4serialization);
  if (sim.dumpVelocityZ) DumpHDF5_MPI2<StreamerVelVectorZ    , Real, FluidGridMPI, LabMPI> (*grid2Dump, sim.time, StreamerVelVectorZ     ::prefix() + name.str(),sim.path4serialization);

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
    #ifdef _DOUBLE_PRECISION_
    ret = ret && 1==fscanf(f, "time: %le\n",   &sim.time);
    ret = ret && 1==fscanf(f, "stepid: %d\n", &sim.step);
    ret = ret && 1==fscanf(f, "uinfx: %le\n", &sim.uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %le\n", &sim.uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %le\n", &sim.uinf[2]);
    #endif
    #ifdef _FLOAT_PRECISION_
    ret = ret && 1==fscanf(f, "time: %e\n",   &sim.time);
    ret = ret && 1==fscanf(f, "stepid: %d\n", &sim.step);
    ret = ret && 1==fscanf(f, "uinfx: %e\n", &sim.uinf[0]);
    ret = ret && 1==fscanf(f, "uinfy: %e\n", &sim.uinf[1]);
    ret = ret && 1==fscanf(f, "uinfz: %e\n", &sim.uinf[2]);
    #endif
    fclose(f);
    if( (not ret) || sim.step<0 || sim.time<0) {
      printf("Error reading restart file. Aborting...\n");
      fflush(0); MPI_Abort(sim.grid->getCartComm(), 1);
    }
  }

  std::stringstream ssR;
  ssR<<"restart_"<<std::setfill('0')<<std::setw(9)<<sim.step;
  if (sim.rank==0) std::cout << "Restarting from " << ssR.str() << "\n";

  ReadHDF5_MPI<StreamerVelocityVector, Real>(* sim.grid, StreamerVelocityVector::prefix()+ssR.str(), sim.path4serialization);

  sim.obstacle_vector->restart(sim.path4serialization+"/"+ssR.str());

  printf("DESERIALIZATION: time is %f and step id is %d\n", sim.time, sim.step);
  // prepare time for next save
  sim.nextSaveTime = sim.time + sim.dumpTime;
}

void Simulation::run()
{
  for (;;) {
    const Real dt = calcMaxTimestep();

    if (timestep(dt)) break;
  }
}

Real Simulation::calcMaxTimestep()
{
  const Real dt_old = sim.dt;
  sim.dt_old = sim.dt;
  const Real hMin = sim.hmin;
  Real CFL = sim.CFL;
  sim.uMax_measured = findMaxU(sim);

  if( CFL > 0 )
  {
    const Real dtDiffusion = (1.0/6.0)*hMin*hMin/(sim.nu+(1.0/6.0)*hMin*sim.uMax_measured);
    const Real dtAdvection = hMin / ( sim.uMax_measured + 1e-8 );
    if ( sim.step < sim.rampup )
    {
      const Real x = sim.step / (Real) sim.rampup;
      const Real rampCFL = std::exp(std::log(1e-3)*(1-x) + std::log(CFL)*x);
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
    printf("=======================================================================\n");
    printf("[CUP3D] step: %d, time: %f, dt: %.2e, uinf: {%f %f %f}, maxU:%f, minH:%f, CFL:%.2e, lambda:%.2e, collision?:%d, blocks:%zu\n",
           sim.step,sim.time, sim.dt, sim.uinf[0], sim.uinf[1], sim.uinf[2],
           sim.uMax_measured, hMin, CFL, sim.lambda, sim.bCollision,
           sim.vInfo().size());
  }

  if (sim.step > sim.step_2nd_start)
  {
    const Real a = dt_old;
    const Real b = sim.dt;
    const Real c1 = -(a+b)/(a*b);
    const Real c2 = b/(a+b)/a;
    sim.coefU[0] = -b*(c1+c2);
    sim.coefU[1] = b*c1;
    sim.coefU[2] = b*c2;
    //sim.coefU[0] = 1.5;
    //sim.coefU[1] = -2.0;
    //sim.coefU[2] = 0.5;
  }
  return sim.dt;
}

bool Simulation::timestep(const Real dt)
{
  const bool bDumpFreq = (sim.saveFreq>0 && (sim.step+ 1)%sim.saveFreq==0);
  const bool bDumpTime = (sim.dumpTime>0 && (sim.time+dt)>sim.nextSaveTime);
  if (bDumpTime) sim.nextSaveTime += sim.dumpTime;
  sim.bDump = (bDumpFreq || bDumpTime);

  //The mesh be adapted before objects are placed on grid
  if (sim.step % 5 == 0 || sim.step < 10) adaptMesh();

  for (size_t c=0; c< sim.pipeline.size() ; c++)
  {
    sim.startProfiler(sim.pipeline[c]->getName());
    (*sim.pipeline[c])(dt);
    sim.stopProfiler();
  }
  sim.step++;
  sim.time += dt;

  if( sim.bDump ) _serialize();

  if (sim.rank == 0 && sim.freqProfiler > 0 && sim.step % sim.freqProfiler == 0)
    sim.printResetProfiler();

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

void Simulation::computeVorticity()
{
  ComputeVorticity findOmega(sim);
  findOmega(0);
}

void Simulation::insertOperator(std::shared_ptr<Operator> op)
{
  sim.pipeline.push_back(std::move(op));
}

void Simulation::touch()
{
  std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < (int)vInfo.size(); ++i)
  {
    const cubism::BlockInfo & info = vInfo[i];
    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    b.clear();
  }
}

CubismUP_3D_NAMESPACE_END
