//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "SimulationData.h"
#include "operators/Operator.h"
#include "obstacles/IF3D_ObstacleVector.h"
#include "Cubism/ArgumentParser.h"

SimulationData::SimulationData(MPI_Comm mpicomm, ArgumentParser &parser) :
  app_comm(mpicomm)
{
  MPI_Comm_rank(app_comm, &rank);
  MPI_Comm_size(app_comm, &nprocs);
  if (rank == 0) parser.print_args();

  // ========== SIMULATION ==========
  // GRID
  bpdx = parser("-bpdx").asInt();
  bpdy = parser("-bpdy").asInt();
  bpdz = parser("-bpdz").asInt();
  nprocsx = parser("-nprocsx").asInt(-1);
  nprocsy = parser("-nprocsy").asInt(-1);
  nprocsz = 1;
  // FLOW
  nu = parser("-nu").asDouble();
  uMax_forced = parser("-uMax_forced").asDouble(0.0);
  lambda = 0.0; //parser("-lambda").asDouble(1e6);
  DLM = 1.0;//parser("-use-dlm").asDouble(0.0);
  CFL = parser("-CFL").asDouble(.1);
  uinf[0] = parser("-uinfx").asDouble(0.0);
  uinf[1] = parser("-uinfy").asDouble(0.0);
  uinf[2] = parser("-uinfz").asDouble(0.0);

  // PIPELINE
  computeDissipation = (bool)parser("-compute-dissipation").asInt(0);

  // OUTPUT
  verbose = parser("-verbose").asBool(true) && rank == 0;
  b2Ddump = parser("-dump2D").asBool(false);
  b3Ddump = parser("-dump3D").asBool(true);

  int dumpFreq = parser("-fdump").asDouble(0);       // dumpFreq==0 means dump freq (in #steps) is not active
  double dumpTime = parser("-tdump").asDouble(0.0);  // dumpTime==0 means dump freq (in time)   is not active
  saveFreq = parser("-fsave").asInt(0);         // dumpFreq==0 means dump freq (in #steps) is not active
  saveTime = parser("-tsave").asDouble(0.0);    // dumpTime==0 means dump freq (in time)   is not active
  rampup = parser("-rampup").asInt(100); // number of dt ramp-up steps

  nsteps = parser("-nsteps").asInt(0);    // 0 to disable this stopping critera.
  endTime = parser("-tend").asDouble(0);  // 0 to disable this stopping critera.

  // TEMP: Removed distinction saving-dumping. Backward compatibility:
  if (saveFreq <= 0 && dumpFreq > 0) saveFreq = dumpFreq;
  if (saveTime <= 0 && dumpTime > 0) saveTime = dumpTime;

  path4serialization = parser("-serialization").asString("./");

  // INITIALIZATION: Mostly unused
  initCond = parser("-initCond").asString("zero");
  useHYPREsolver = parser("-useHYPREsolver").asString("");
  // BOUNDARY CONDITIONS
  // accepted dirichlet, periodic, freespace/unbounded, fakeOpen
  std::string BC_x = parser("-BC_x").asString("dirichlet");
  std::string BC_y = parser("-BC_y").asString("dirichlet");
  std::string BC_z = parser("-BC_z").asString("dirichlet");
  const Real fadeLen = parser("-fade_len").asDouble(.005);
  // BC
  if(BC_x=="unbounded") BC_x = "freespace"; // tomato tomato
  if(BC_y=="unbounded") BC_y = "freespace"; // tomato tomato
  if(BC_z=="unbounded") BC_z = "freespace"; // tomato tomato
  // boundary killing useless for unbounded or periodic
  fadeOutLengthPRHS[0] = BC_x=="dirichlet"? fadeLen : 0;
  fadeOutLengthPRHS[1] = BC_y=="dirichlet"? fadeLen : 0;
  fadeOutLengthPRHS[2] = BC_z=="dirichlet"? fadeLen : 0;

  if(BC_x=="freespace" || BC_y=="freespace" || BC_z=="freespace")
  {
    if(BC_x=="freespace" && BC_y=="freespace" && BC_z=="freespace") {
      bUseUnboundedBC = true; // poisson solver
    } else {
     fprintf(stderr,"ERROR: either all or no BC can be freespace/unbounded!\n");
     fflush(0); abort();
    }
  }

  if(BC_x=="fakeopen" || BC_y=="fakeopen" || BC_z=="fakeopen")
  {
    if(BC_x=="fakeopen" && BC_y=="fakeopen" && BC_z=="fakeopen")
    {
      fadeOutLengthU[0] = fadeLen; fadeOutLengthPRHS[0] = fadeLen;
      fadeOutLengthU[1] = fadeLen; fadeOutLengthPRHS[1] = fadeLen;
      fadeOutLengthU[2] = fadeLen; fadeOutLengthPRHS[2] = fadeLen;
      BC_x = "freespace"; BC_y = "freespace"; BC_z = "freespace";
      bUseFourierBC = true; // poisson solver
    } else {
     fprintf(stderr,"ERROR: either all or no BC can be fakeopen!\n");
     fflush(0); abort();
    }
  }

  BCx_flag = string2BCflag(BC_x);
  BCy_flag = string2BCflag(BC_y);
  BCz_flag = string2BCflag(BC_z);


  // DFT if we are periodic in all directions:
  if(BC_x=="periodic"&&BC_y=="periodic"&&BC_z=="periodic") bUseFourierBC = true;
  if(rank==0)
  printf("Boundary pressure RHS / FD smoothing region sizes {%f,%f,%f}\n",
    fadeOutLengthPRHS[0], fadeOutLengthPRHS[1], fadeOutLengthPRHS[2]);
  _argumentsSanityCheck();

  // ============ REST =============
}

SimulationData::SimulationData(MPI_Comm mpicomm) : app_comm(mpicomm) {}

void SimulationData::_argumentsSanityCheck()
{
  // Grid.
  if (bpdx < 1 || bpdy < 1 || bpdz < 1) {
      fprintf(stderr, "Invalid bpd: %d x %d x %d\n", bpdx, bpdy, bpdz);
      abort();
  }

  // Flow.
  assert(nu >= 0);
  assert(lambda > 0 || DLM >= 1);
  assert(CFL > 0.0);

  // Output.
  assert(saveFreq >= 0.0);
  assert(saveTime >= 0.0);
}

SimulationData::~SimulationData()
{
  delete grid;
  delete obstacle_vector;
  while(!pipeline.empty()) {
    auto * g = pipeline.back();
    pipeline.pop_back();
    delete g;
  }
  #ifdef CUP_ASYNC_DUMP
    if(dumper not_eq nullptr) {
      dumper->join();
      delete dumper;
    }
    delete dump;
    MPI_Comm_free(&dump_comm);
  #endif
}

void SimulationData::startProfiler(std::string name) const
{
  #ifndef SMARTIES_APP
    profiler->push_start(name);
  #endif
}
void SimulationData::stopProfiler() const
{
  #ifndef SMARTIES_APP
    profiler->pop_stop();
  #endif
}
void SimulationData::printResetProfiler()
{
  #ifndef SMARTIES_APP
    profiler->printSummary();
    profiler->reset();
  #endif
}
