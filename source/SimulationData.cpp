//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include <unistd.h>

#include "SimulationData.h"
#include "operators/Operator.h"
#include "Obstacles/ObstacleVector.h"
#include <Cubism/ArgumentParser.h>
#include <Cubism/Profiler.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;
  
BCflag cubismBCX; 
BCflag cubismBCY; 
BCflag cubismBCZ; 
SimulationData::SimulationData(MPI_Comm mpicomm, ArgumentParser &parser): comm(mpicomm)
{
  // Initialize MPI related variables
  MPI_Comm_rank(comm, &rank);

  // Print parser content
  if (rank == 0) parser.print_args();

  // ========== PARSE ARGUMENTS ==========

  // restart the simulation?
  bRestart = parser("-restart").asBool(false);
  checkpoint_steps = parser("-checkpointsteps").asInt(1000);

  // BLOCKS PER DIMENSION
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();

  // AMR SETTINGS
  levelMax = parser("-levelMax").asInt();
  levelStart = parser("-levelStart").asInt(levelMax-1);
  Rtol = parser("-Rtol").asDouble();
  Ctol = parser("-Ctol").asDouble();
  levelMaxVorticity = parser("-levelMaxVorticity").asInt(levelMax);

  // SIMULATION DOMAIN
  extents[0] = parser("extentx").asDouble(0);
  extents[1] = parser("extenty").asDouble(0);
  extents[2] = parser("extentz").asDouble(0);
  if (extents[0] + extents[1] + extents[2] < 1e-21) extents[0] = parser("extent").asDouble(1);

  // SPEED OF FRAME OF REFERENCE
  uinf[0] = parser("-uinfx").asDouble(0.0);
  uinf[1] = parser("-uinfy").asDouble(0.0);
  uinf[2] = parser("-uinfz").asDouble(0.0);

  // TIMESTEPPING
  CFL = parser("-CFL").asDouble(.1);
  dt = parser("-dt").asDouble(0);
  rampup = parser("-rampup").asInt(100); // number of dt ramp-up steps
  nsteps = parser("-nsteps").asInt(0);    // 0 to disable this stopping critera.
  endTime = parser("-tend").asDouble(0);  // 0 to disable this stopping critera.
  step_2nd_start = 2;

  // FLOW
  nu = parser("-nu").asDouble();

  // IC
  initCond = parser("-initCond").asString("zero");

  // SPEED FOR CHANNEL FLOW
  uMax_forced = parser("-uMax_forced").asDouble(0.0);
  bFixMassFlux = parser("-bFixMassFlux").asBool(false);

  // PENALIZATION
  bImplicitPenalization = parser("-implicitPenalization").asBool(true);
  lambda = parser("-lambda").asDouble(1e6);
  DLM = parser("-use-dlm").asDouble(0);

  // DISSIPATION DIAGNOSTIC
  freqDiagnostics = parser("-freqDiagnostics").asInt(100);

  // PROFILER
  freqProfiler = parser("-freqProfiler").asInt(0);

  // POISSON SOLVER
  PoissonErrorTol = parser("-poissonTol").asDouble(1e-6); // absolute error
  PoissonErrorTolRel = parser("-poissonTolRel").asDouble(1e-4); // relative error
  bMeanConstraint = parser("-bMeanConstraint").asInt(1); //zero mean constraint 
  poissonSolver = parser("-poissonSolver").asString("iterative");

  // IMPLICIT DIFFUSION SOLVER
  implicitDiffusion = parser("-implicitDiffusion").asBool(false);
  DiffusionErrorTol = parser("-diffusionTol").asDouble(1e-6); // absolute error
  DiffusionErrorTolRel = parser("diffusionTolRel").asDouble(1e-4); // relative error

  uMax_allowed = parser("-umax").asDouble(10.0);

  // BOUNDARY CONDITIONS
  // accepted periodic, freespace or wall
  std::string BC_x = parser("-BC_x").asString("freespace");
  std::string BC_y = parser("-BC_y").asString("freespace");
  std::string BC_z = parser("-BC_z").asString("freespace");
  BCx_flag = string2BCflag(BC_x);
  BCy_flag = string2BCflag(BC_y);
  BCz_flag = string2BCflag(BC_z);

  cubismBCX = BCx_flag;
  cubismBCY = BCy_flag;
  cubismBCZ = BCz_flag;

  // OUTPUT
  muteAll = parser("-muteAll").asInt(0);
  verbose = muteAll ? false : parser("-verbose").asInt(1) && rank == 0;
  int dumpFreq = parser("-fdump").asDouble(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  dumpTime = parser("-tdump").asDouble(0.0);    // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);         // dumpFreq==0 means dump freq (in #steps) is not active

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if (saveFreq <= 0 && dumpFreq > 0) saveFreq = dumpFreq;
  path4serialization = parser("-serialization").asString("./");

  // Dumping
  dumpChi       = parser("-dumpChi"      ).asBool(true);
  dumpOmega     = parser("-dumpOmega"    ).asBool(true);
  dumpP         = parser("-dumpP"        ).asBool(false);
  dumpOmegaX    = parser("-dumpOmegaX"   ).asBool(false);
  dumpOmegaY    = parser("-dumpOmegaY"   ).asBool(false);
  dumpOmegaZ    = parser("-dumpOmegaZ"   ).asBool(false);
  dumpVelocity  = parser("-dumpVelocity" ).asBool(false);
  dumpVelocityX = parser("-dumpVelocityX").asBool(false);
  dumpVelocityY = parser("-dumpVelocityY").asBool(false);
  dumpVelocityZ = parser("-dumpVelocityZ").asBool(false);
}

void SimulationData::_preprocessArguments()
{
  assert(profiler == nullptr);  // This should not be possible at all.
  profiler = new cubism::Profiler();
  if (bpdx < 1 || bpdy < 1 || bpdz < 1)
  {
    fprintf(stderr, "Invalid bpd: %d x %d x %d\n", bpdx, bpdy, bpdz);
    fflush(0); abort();
  }
  const int aux = 1 << (levelMax -1);
  const Real NFE[3] = {
      (Real) bpdx * aux * ScalarBlock::sizeX,
      (Real) bpdy * aux * ScalarBlock::sizeY,
      (Real) bpdz * aux * ScalarBlock::sizeZ,
  };
  const Real maxbpd = std::max({NFE[0], NFE[1], NFE[2]});
  maxextent = std::max({extents[0], extents[1], extents[2]});
  if( extents[0] <= 0 || extents[1] <= 0 || extents[2] <= 0 )
  {
    extents[0] = (NFE[0]/maxbpd) * maxextent;
    extents[1] = (NFE[1]/maxbpd) * maxextent;
    extents[2] = (NFE[2]/maxbpd) * maxextent;
  }
  else
  {
    fprintf(stderr, "Invalid extent: %f x %f x %f\n", extents[0], extents[1], extents[2]);
    fflush(0); abort();
  }
  hmin = extents[0] / NFE[0];
  hmax = extents[0] * aux / NFE[0];
  assert(nu >= 0);
  assert(lambda > 0 || DLM > 0);
  assert(saveFreq >= 0.0);
  assert(dumpTime >= 0.0);
}

SimulationData::~SimulationData()
{
  delete profiler;
  delete obstacle_vector;
  delete chi;
  delete vel;
  delete lhs;
  delete tmpV;
  delete pres;
  delete chi_amr;
  delete vel_amr;
  delete lhs_amr;
  delete tmpV_amr;
  delete pres_amr;
}

void SimulationData::startProfiler(std::string name) const
{
  profiler->push_start(name);
}
void SimulationData::stopProfiler() const
{
  profiler->pop_stop();
}
void SimulationData::printResetProfiler()
{
  profiler->printSummary();
  profiler->reset();
}

void SimulationData::writeRestartFiles()
{
  // write restart file for field
  if (rank == 0)
  {
     std::stringstream ssR;
     ssR << path4serialization + "/field.restart";
     FILE * fField = fopen(ssR.str().c_str(), "w");
     if (fField == NULL)
     {
        printf("Could not write %s. Aborting...\n", "field.restart");
        fflush(0); abort();
     }
     assert(fField != NULL);
     fprintf(fField, "time: %20.20e\n",  (double)time);
     fprintf(fField, "stepid: %d\n",     step);
     fprintf(fField, "uinfx: %20.20e\n", (double)uinf[0]);
     fprintf(fField, "uinfy: %20.20e\n", (double)uinf[1]);
     fprintf(fField, "uinfz: %20.20e\n", (double)uinf[2]);
     fprintf(fField, "dt: %20.20e\n",    (double)dt);
     fclose(fField);
  }

  // write restart file for shapes
  int size;
  MPI_Comm_size(comm,&size);
  const size_t tasks = obstacle_vector->nObstacles();
  size_t my_share = tasks / size;
  if (tasks % size != 0 && rank == size - 1) //last rank gets what's left
  {
   my_share += tasks % size;
  }
  const size_t my_start = rank * (tasks/ size);
  const size_t my_end   = my_start + my_share;

  #pragma omp parallel for schedule(static,1)
  for(size_t j = my_start ; j < my_end ; j++)
  {
    auto & shape = obstacle_vector->getObstacleVector()[j];
    std::stringstream ssR;
    ssR << path4serialization + "/shape_" << shape->obstacleID << ".restart";
    FILE * fShape = fopen(ssR.str().c_str(), "w");
    if (fShape == NULL)
    {
      printf("Could not write %s. Aborting...\n", ssR.str().c_str());
      fflush(0); abort();
    }
    shape->saveRestart( fShape );
    fclose(fShape);
  }
}

void SimulationData::readRestartFiles()
{
  // read restart file for field
  FILE * fField = fopen("field.restart", "r");
  if (fField == NULL) {
    printf("Could not read %s. Aborting...\n", "field.restart");
    fflush(0); abort();
  }
  assert(fField != NULL);
  if (rank == 0 && verbose) printf("Reading %s...\n", "field.restart");
  bool ret = true;
  double in_time, in_uinfx, in_uinfy, in_uinfz, in_dt;
  ret = ret && 1==fscanf(fField, "time: %le\n",   &in_time);
  ret = ret && 1==fscanf(fField, "stepid: %d\n",  &step);
  ret = ret && 1==fscanf(fField, "uinfx: %le\n",  &in_uinfx);
  ret = ret && 1==fscanf(fField, "uinfy: %le\n",  &in_uinfy);
  ret = ret && 1==fscanf(fField, "uinfz: %le\n",  &in_uinfz);
  ret = ret && 1==fscanf(fField, "dt: %le\n",     &in_dt);
  time  = (Real) in_time ;
  uinf[0] = (Real) in_uinfx;
  uinf[1] = (Real) in_uinfy;
  uinf[2] = (Real) in_uinfz;
  dt    = (Real) in_dt   ;
  fclose(fField);
  if( (not ret) || step<0 || time<0) {
    printf("Error reading restart file. Aborting...\n");
    fflush(0); abort();
  }
  if (rank == 0 && verbose) printf("Restarting flow.. time: %le, stepid: %d, uinfx: %le, uinfy: %le, uinfz: %le\n", (double)time, step, (double)uinf[0], (double)uinf[1], (double)uinf[2]);
  nextSaveTime = time + dumpTime;

  // read restart file for shapes
  for (auto &shape : obstacle_vector->getObstacleVector())
  {
    std::stringstream ssR;
    ssR << "shape_" << shape->obstacleID << ".restart";
    FILE * fShape = fopen(ssR.str().c_str(), "r");
    if (fShape == NULL) {
      printf("Could not read %s. Aborting...\n", ssR.str().c_str());
      fflush(0); abort();
    }
    if (rank == 0 && verbose) printf("Reading %s...\n", ssR.str().c_str());
    shape->loadRestart( fShape );
    fclose(fShape);
  }
}

CubismUP_3D_NAMESPACE_END
