//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "Simulation.h"

#include "operators/InitialConditions.h"
#include "operators/ObstaclesCreate.h"
#include "operators/AdvectionDiffusion.h"
#include "operators/ObstaclesUpdate.h"
#include "operators/Penalization.h"
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

Simulation::Simulation(MPI_Comm mpicomm, ArgumentParser & parser) : sim(mpicomm, parser)
{
  if( sim.verbose )
  {
    #pragma omp parallel
    {
          int numThreads = omp_get_num_threads();
	  int size;
	  MPI_Comm_size(mpicomm,&size);
          #pragma omp master
	  std::cout << "[CUP3D] Running with " << size  << " rank(s) and " << numThreads <<  " thread(s)." << std::endl;
    }
  }
  // Make sure given arguments are valid
  if( sim.verbose )
    std::cout << "[CUP3D] Parsing Arguments.. " << std::endl;
  sim._preprocessArguments();

  // Setup and Initialize Grid
  if( sim.verbose )
    std::cout << "[CUP3D] Allocating Grid.. " << std::endl;
  setupGrid();

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
    deserialize();
  else
  {
    _ic();
    if( sim.verbose )
      std::cout << "[CUP3D] Performing Initial Refinement of Grid.. " << std::endl;
    initialGridRefinement();
  }
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
  }
}

void Simulation::adaptMesh()
{
  sim.startProfiler("Mesh refinement");

  computeVorticity();
  compute<ScalarLab>(GradChiOnTmp(sim),sim.chi);

  sim.tmpV_amr->Tag();
  sim.lhs_amr ->TagLike(sim.tmpVInfo());
  sim.vel_amr ->TagLike(sim.tmpVInfo());
  sim.chi_amr ->TagLike(sim.tmpVInfo());
  sim.pres_amr->TagLike(sim.tmpVInfo());
  sim.chi_amr ->Adapt(sim.time,sim.verbose,true);
  sim.lhs_amr ->Adapt(sim.time,false,true);
  sim.tmpV_amr->Adapt(sim.time,false,true);
  sim.pres_amr->Adapt(sim.time,false,false);
  sim. vel_amr->Adapt(sim.time,false,false);

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

void Simulation::setupGrid()
{
  sim.chi  = new ScalarGrid(sim.bpdx,sim.bpdy,sim.bpdz,sim.maxextent,sim.levelStart,sim.levelMax,sim.comm,(sim.BCx_flag == periodic),(sim.BCy_flag == periodic),(sim.BCz_flag == periodic));
  sim.lhs  = new ScalarGrid(sim.bpdx,sim.bpdy,sim.bpdz,sim.maxextent,sim.levelStart,sim.levelMax,sim.comm,(sim.BCx_flag == periodic),(sim.BCy_flag == periodic),(sim.BCz_flag == periodic));
  sim.pres = new ScalarGrid(sim.bpdx,sim.bpdy,sim.bpdz,sim.maxextent,sim.levelStart,sim.levelMax,sim.comm,(sim.BCx_flag == periodic),(sim.BCy_flag == periodic),(sim.BCz_flag == periodic));
  sim.vel  = new VectorGrid(sim.bpdx,sim.bpdy,sim.bpdz,sim.maxextent,sim.levelStart,sim.levelMax,sim.comm,(sim.BCx_flag == periodic),(sim.BCy_flag == periodic),(sim.BCz_flag == periodic));
  sim.tmpV = new VectorGrid(sim.bpdx,sim.bpdy,sim.bpdz,sim.maxextent,sim.levelStart,sim.levelMax,sim.comm,(sim.BCx_flag == periodic),(sim.BCy_flag == periodic),(sim.BCz_flag == periodic));
  //Refine/compress only according to chi field for now
  sim.chi_amr  = new ScalarAMR( *(sim.chi ),sim.Rtol,sim.Ctol);
  sim.lhs_amr  = new ScalarAMR( *(sim.lhs ),sim.Rtol,sim.Ctol);
  sim.pres_amr = new ScalarAMR( *(sim.pres),sim.Rtol,sim.Ctol);
  sim.vel_amr  = new VectorAMR( *(sim.vel ),sim.Rtol,sim.Ctol);
  sim.tmpV_amr = new VectorAMR( *(sim.tmpV),sim.Rtol,sim.Ctol);
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
  // overwrites tmpU, tmpV, tmpW and pressure solver's RHS. Then,
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

void Simulation::serialize(const std::string append)
{
  sim.startProfiler("DumpHDF5_MPI");

  std::stringstream name;
  if (append == "") name<<"_";
  else name<<append;
  name<<std::setfill('0')<<std::setw(9)<<sim.step;

  if (sim.dumpOmega || sim.dumpOmegaX || sim.dumpOmegaY || sim.dumpOmegaZ)
    computeVorticity();

  //dump multi-block datasets with scalar quantities or magnitude of vector quantities
  if (sim.dumpP        ) DumpHDF5_MPI2<cubism::StreamerScalar,Real,ScalarGrid> (*sim.pres, sim.time, "pres" + name.str(),sim.path4serialization);
  if (sim.dumpChi      ) DumpHDF5_MPI2<cubism::StreamerScalar,Real,ScalarGrid> (*sim.chi , sim.time, "chi"  + name.str(),sim.path4serialization);
  if (sim.dumpOmega    ) DumpHDF5_MPI2<cubism::StreamerVector,Real,VectorGrid> (*sim.tmpV, sim.time, "tmp"  + name.str(),sim.path4serialization);
  if (sim.dumpVelocity ) DumpHDF5_MPI2<cubism::StreamerVector,Real,VectorGrid> (*sim.vel , sim.time, "vel"  + name.str(),sim.path4serialization);

  //dump components of vectors
  if (sim.dumpOmegaX   ) DumpHDF5_MPI2<StreamerVectorX,Real,VectorGrid> (*sim.tmpV, sim.time, "tmpX" + name.str(),sim.path4serialization);
  if (sim.dumpOmegaY   ) DumpHDF5_MPI2<StreamerVectorY,Real,VectorGrid> (*sim.tmpV, sim.time, "tmpY" + name.str(),sim.path4serialization);
  if (sim.dumpOmegaZ   ) DumpHDF5_MPI2<StreamerVectorZ,Real,VectorGrid> (*sim.tmpV, sim.time, "tmpZ" + name.str(),sim.path4serialization);
  if (sim.dumpVelocityX) DumpHDF5_MPI2<StreamerVectorX,Real,VectorGrid> (*sim.vel, sim.time, "velX" + name.str(),sim.path4serialization);
  if (sim.dumpVelocityY) DumpHDF5_MPI2<StreamerVectorY,Real,VectorGrid> (*sim.vel, sim.time, "velY" + name.str(),sim.path4serialization);
  if (sim.dumpVelocityZ) DumpHDF5_MPI2<StreamerVectorZ,Real,VectorGrid> (*sim.vel, sim.time, "velZ" + name.str(),sim.path4serialization);

  sim.stopProfiler();
}

void Simulation::deserialize()
{
  // create filename from step
  sim.readRestartFiles();
  std::stringstream ss;
  ss<<"restart_"<<std::setfill('0')<<std::setw(9)<<sim.step;

  const std::vector<BlockInfo>&  chiInfo = sim.chi ->getBlocksInfo();
  const std::vector<BlockInfo>&  velInfo = sim.vel ->getBlocksInfo();
  const std::vector<BlockInfo>& tmpVInfo = sim.tmpV->getBlocksInfo();
  const std::vector<BlockInfo>&  lhsInfo = sim.lhs ->getBlocksInfo();

  //The only field that is needed for restarting is velocity. Chi is derived from the files we
  //read for obstacles. Here we also read pres so that the Poisson solver has the same
  //initial guess, which in turn leads to restarted simulations having the exact same result
  //as non-restarted ones (we also read pres because we need to read at least
  //one ScalarGrid, see hack below).
  ReadHDF5_MPI<StreamerVector, Real>(*(sim.vel ), "vel_"  + ss.str(), sim.path4serialization);
  ReadHDF5_MPI<StreamerScalar, Real>(*(sim.pres), "pres_" + ss.str(), sim.path4serialization);

  //hack: need to "read" the other grids too, so that the mesh is the same for every grid.
  //So we read VectorGrids from "vel" and ScalarGrids from "pres". We don't care about the
  //grid point values (those are set to zero below), we only care about the grid structure,
  //i.e. refinement levels etc.
  ReadHDF5_MPI<StreamerScalar, Real>(*(sim.chi ), "pres_" + ss.str(), sim.path4serialization);
  ReadHDF5_MPI<StreamerScalar, Real>(*(sim.lhs ), "pres_" + ss.str(), sim.path4serialization);
  ReadHDF5_MPI<StreamerVector, Real>(*(sim.tmpV),  "vel_" + ss.str(), sim.path4serialization);
  #pragma omp parallel for
  for (size_t i=0; i < velInfo.size(); i++)
  {
    ScalarBlock& CHI  = *(ScalarBlock*)  chiInfo[i].ptrBlock;  CHI.clear();
    ScalarBlock& LHS  = *(ScalarBlock*)  lhsInfo[i].ptrBlock;  LHS.clear();
    VectorBlock& TMPV = *(VectorBlock*) tmpVInfo[i].ptrBlock; TMPV.clear();
  }
}

void Simulation::simulate()
{
  for (;;) {
    const Real dt = calcMaxTimestep();

    if (advance(dt)) break;
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
    fflush(0); MPI_Abort(sim.comm, 1);
  }


  // if DLM>0, adapt lambda such that penal term is independent of time step
  if (sim.DLM > 0) sim.lambda = sim.DLM / sim.dt;

  if( sim.rank == 0 ) {
    printf("=======================================================================\n");
    printf("[CUP3D] step: %d, time: %f, dt: %.2e, uinf: {%f %f %f}, maxU:%f, minH:%f, CFL:%.2e, lambda:%.2e, collision?:%d, blocks:%zu\n",
           sim.step,sim.time, sim.dt, sim.uinf[0], sim.uinf[1], sim.uinf[2],
           sim.uMax_measured, hMin, CFL, sim.lambda, sim.bCollision,
           sim.velInfo().size());
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

bool Simulation::advance(const Real dt)
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

  if( sim.bDump ) serialize();
  if (sim.step % sim.checkpoint_steps == 0)  //checkpoint for restarting
  {
    std::stringstream name;
    name<<"restart_"<<std::setfill('0')<<std::setw(9)<<sim.step;
    DumpHDF5_MPI<StreamerScalar, Real> (*sim.pres, sim.time, "pres_" + name.str(),sim.path4serialization, false);
    DumpHDF5_MPI<StreamerVector, Real> (*sim.vel , sim.time, "vel_"  + name.str(),sim.path4serialization, false);
    sim.writeRestartFiles();
  }

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

CubismUP_3D_NAMESPACE_END
