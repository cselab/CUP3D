//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include <unistd.h>

#include "SimulationData.h"
#include "operators/Operator.h"
#include "obstacles/ObstacleVector.h"
#include <Cubism/ArgumentParser.h>
#include <Cubism/Profiler.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

SimulationData::SimulationData(MPI_Comm mpicomm, ArgumentParser &parser): app_comm(mpicomm)
{
  // Initialize MPI related variables
  MPI_Comm_rank(app_comm, &rank);
  MPI_Comm_size(app_comm, &nprocs);

  // Print parser content
  if (rank == 0) parser.print_args();

  // ========== PARSE ARGUMENTS ==========
  // BLOCKS PER DIMENSION
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();

  // AMR SETTINGS
  levelMax = parser("-levelMax").asInt();
  levelStart = parser("-levelStart").asInt(levelMax-1);
  Rtol = parser("-Rtol").asDouble();
  Ctol = parser("-Ctol").asDouble();

  // SIMULATION DOMAIN
  extent[0] = parser("extentx").asDouble(1);
  extent[1] = parser("extenty").asDouble(0);
  extent[2] = parser("extentz").asDouble(0);

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
  icFromH5 = parser("-icFromH5").asString("");

  // SPEED FOR CHANNEL FLOW
  uMax_forced = parser("-uMax_forced").asDouble(0.0);

  // PENALIZATION
  bImplicitPenalization = parser("-implicitPenalization").asBool(true);
  lambda = parser("-lambda").asDouble(1e6);
  DLM = parser("-use-dlm").asDouble(0);

  // DISSIPATION DIAGNOSTIC
  freqDiagnostics = parser("-freqDiagnostics").asInt(100);

  // POISSON SOLVER
  PoissonErrorTol = parser("-poissonTol").asDouble(1e-6); // absolute error
  PoissonErrorTolRel = parser("-poissonTolRel").asDouble(1e-4); // relative error

  // BOUNDARY CONDITIONS
  // accepted periodic, freespace or wall
  std::string BC_x = parser("-BC_x").asString("freespace");
  std::string BC_y = parser("-BC_y").asString("freespace");
  std::string BC_z = parser("-BC_z").asString("freespace");
  BCx_flag = string2BCflag(BC_x);
  BCy_flag = string2BCflag(BC_y);
  BCz_flag = string2BCflag(BC_z);

  // OUTPUT
  muteAll = parser("-muteAll").asInt(0);
  verbose = muteAll ? false : parser("-verbose").asInt(1) && rank == 0;
  int dumpFreq = parser("-fdump").asDouble(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  double dumpTime = parser("-tdump").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);         // dumpFreq==0 means dump freq (in #steps) is not active
  saveTime = parser("-tsave").asDouble(0.0);    // dumpTime==0 means dump freq (in time)   is not active

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if (saveFreq <= 0 && dumpFreq > 0) saveFreq = dumpFreq;
  if (saveTime <= 0 && dumpTime > 0) saveTime = dumpTime;
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
  int aux = 1 << (levelMax -1);
  const double NFE[3] = {
      (double) bpdx * aux * FluidBlock::sizeX,
      (double) bpdy * aux * FluidBlock::sizeY,
      (double) bpdz * aux * FluidBlock::sizeZ,
  };
  const double maxbpd = std::max({NFE[0], NFE[1], NFE[2]});
  maxextent = std::max({extent[0], extent[1], extent[2]});
  if( extent[0] <= 0 || extent[1] <= 0 || extent[2] <= 0 )
  {
    extent[0] = (NFE[0]/maxbpd) * maxextent;
    extent[1] = (NFE[1]/maxbpd) * maxextent;
    extent[2] = (NFE[2]/maxbpd) * maxextent;
  }
  else
  {
    fprintf(stderr, "Invalid extent: %f x %f x %f\n", extent[0], extent[1], extent[2]);
    fflush(0); abort();
  }
  hmin = extent[0] / NFE[0];
  hmax = extent[0] * aux / NFE[0];
  assert(nu >= 0);
  assert(lambda > 0 || DLM > 0);
  assert(saveFreq >= 0.0);
  assert(saveTime >= 0.0);
}

SimulationData::~SimulationData()
{
  delete grid;
  delete gridPoisson;
  delete profiler;
  delete obstacle_vector;
  delete amr;
  delete amr2;
  while(!pipeline.empty()) {
    auto * g = pipeline.back();
    pipeline.pop_back();
    delete g;
  }
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

CubismUP_3D_NAMESPACE_END
