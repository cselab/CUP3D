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
  MPI_Comm_rank(app_comm, &rank);
  MPI_Comm_size(app_comm, &nprocs);

  if (rank == 0) parser.print_args();

  // ========== SIMULATION ==========
  // GRID
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();
  levelMax = parser("-levelMax").asInt(1);
  levelStart = parser("-levelStart").asInt(levelMax-1);
  Rtol = parser("-Rtol").asDouble(1.0);
  Ctol = parser("-Ctol").asDouble(0.1);

  extent[0] = parser("extentx").asDouble(1);
  extent[1] = parser("extenty").asDouble(0);
  extent[2] = parser("extentz").asDouble(0);

  // FLOW
  nu = parser("-nu").asDouble();

  // IC
  initCond = parser("-initCond").asString("zero");
  spectralIC = parser("-spectralIC").asString("");
  k0 = parser("-k0").asDouble(10.0);
  tke0 = parser("-tke0").asDouble(1.0);
  icFromH5 = parser("-icFromH5").asString("");

  // FORCING HIT=
  turbKinEn_target = parser("-turbKinEn_target").asDouble(0);
  enInjectionRate = parser("-energyInjectionRate").asDouble(0);
  const bool bSpectralForcingHint = turbKinEn_target>0 || enInjectionRate>0;
  spectralForcing = parser("-spectralForcing").asBool(bSpectralForcingHint);
  if(turbKinEn_target>0 && enInjectionRate>0) {
    fprintf(stderr,"ERROR: either constant energy injection rate "
                   "or forcing to fixed energy target\n");
    fflush(0); abort();
  }

  // PIPELINE && FORCING
  freqDiagnostics = parser("-freqDiagnostics").asInt(100);
  bImplicitPenalization = parser("-implicitPenalization").asBool(false);
  bChannelFixedMassFlux = parser("-channelFixedMassFlux").asBool(false);

  uMax_forced = parser("-uMax_forced").asDouble(0.0);

  // SGS
  sgs = parser("-sgs").asString("");
  cs = parser("-cs").asDouble(0.2);
  bComputeCs2Spectrum = parser("-cs2spectrum").asBool(false);

  // SGS_RL
  sgs_rl = parser("-sgs_rl").asBool(false);
  nAgentsPerBlock = parser("-nAgentsPerBlock").asInt(1);

  lambda = parser("-lambda").asDouble(1e6);
  DLM = parser("-use-dlm").asDouble(0);
  CFL = parser("-CFL").asDouble(.1);
  dt = parser("-dt").asDouble(0);
  uinf[0] = parser("-uinfx").asDouble(0.0);
  uinf[1] = parser("-uinfy").asDouble(0.0);
  uinf[2] = parser("-uinfz").asDouble(0.0);

  // OUTPUT
  verbose = parser("-verbose").asBool(true) && rank == 0;

  // ANALYSIS
  analysis = parser("-analysis").asString("");
  timeAnalysis = parser("-tAnalysis").asDouble(0.0);
  freqAnalysis = parser("-fAnalysis").asInt(0);

  int dumpFreq = parser("-fdump").asDouble(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  double dumpTime = parser("-tdump").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);         // dumpFreq==0 means dump freq (in #steps) is not active
  saveTime = parser("-tsave").asDouble(0.0);    // dumpTime==0 means dump freq (in time)   is not active
  rampup = parser("-rampup").asInt(100); // number of dt ramp-up steps
  PoissonErrorTol = parser("-poissonTol").asDouble(1e-6); // absolute error tolerance for poisson eq.
  PoissonErrorTolRel = parser("-poissonTolRel").asDouble(1e-4); // relative error tolerance for poisson eq.

  nsteps = parser("-nsteps").asInt(0);    // 0 to disable this stopping critera.
  endTime = parser("-tend").asDouble(0);  // 0 to disable this stopping critera.

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if (saveFreq <= 0 && dumpFreq > 0) saveFreq = dumpFreq;
  if (saveTime <= 0 && dumpTime > 0) saveTime = dumpTime;

  path4serialization = parser("-serialization").asString("./");

  // INITIALIZATION: Mostly unused
  useSolver = parser("-useSolver").asString("");
  // BOUNDARY CONDITIONS
  // accepted dirichlet, periodic, freespace/unbounded, fakeOpen
  std::string BC_x = parser("-BC_x").asString("dirichlet");
  std::string BC_y = parser("-BC_y").asString("dirichlet");
  std::string BC_z = parser("-BC_z").asString("dirichlet");
  BCx_flag = string2BCflag(BC_x);
  BCy_flag = string2BCflag(BC_y);
  BCz_flag = string2BCflag(BC_z);

  // order of accuracy of timestepping
  TimeOrder = parser("-TimeOrder").asInt(1);
  assert (TimeOrder == 1 || TimeOrder == 2);
  step_2nd_start = 50;

  //Dumping
  dumpChi = parser("-dumpChi").asBool(true);
  dumpOmega = parser("-dumpOmega").asBool(true);
  dumpP = parser("-dumpP").asBool(false);
  dumpOmegaX = parser("-dumpOmegaX").asBool(false);
  dumpOmegaY = parser("-dumpOmegaY").asBool(false);
  dumpOmegaZ = parser("-dumpOmegaZ").asBool(false);
  dumpVelocity = parser("-dumpVelocity").asBool(false);
  dumpVelocityX = parser("-dumpVelocityX").asBool(false);
  dumpVelocityY = parser("-dumpVelocityY").asBool(false);
  dumpVelocityZ = parser("-dumpVelocityZ").asBool(false);

  // ============ REST =============
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
    abort();
  }
  assert(nu >= 0);
  assert(lambda > 0 || DLM > 0);
  assert(saveFreq >= 0.0);
  assert(saveTime >= 0.0);
}

SimulationData::~SimulationData()
{
  delete grid;
  delete profiler;
  delete obstacle_vector;

  if(amr not_eq nullptr) delete amr;

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
