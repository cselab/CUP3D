//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#pragma once

#include "Definitions.h"
#include "Profiler.h"

class Shape;

struct SimulationData
{
  #ifndef SMARTIES_APP
    Profiler * profiler = new Profiler();
  #endif
  VectorGrid * U  = nullptr;
  VectorGrid * Us = nullptr;
  VectorGrid * F  = nullptr;
  VectorGrid * V  = nullptr;

  ScalarGrid * pRHS = nullptr;
  ScalarGrid * chi  = nullptr;
  ScalarGrid * pres = nullptr;
  ScalarGrid * tmp  = nullptr;

  const inline std::vector<BlockInfo>& uInfo() {
    return U->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& sInfo() {
    return Us->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& fInfo() {
    return F->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& vInfo() {
    return V->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& rInfo() {
    return pRHS->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& cInfo() {
    return chi->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& pInfo() {
    return pres->getBlocksInfo();
  }
  const inline std::vector<BlockInfo>& tInfo() {
    return tmp->getBlocksInfo();
  }

  void allocateGrid();

  std::vector<Shape*> shapes;
  // simulation status
  // nsteps==0 means that this stopping criteria is not active
  int step=0, nsteps=0;
  // endTime==0  means that this stopping criteria is not active
  double time=0, endTime=0;
  double dt = 0;
  // mpi
  int rank=-1, nprocs=-1;
  // grid
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;
  // flow variables
  Real uinf[3] = {0, 0, 0};
  double nu=0, CFL=0, lambda=-1, DLM=1;
  // simulation settings
  bool computeDissipation=false;
  bool b3Ddump=true, b2Ddump=false, bDump=false;
  bool rampup=true;
  bool verbose=false;
  bool muteAll = false;
  #ifndef CUP_UNBOUNDED_FFT
    Real fadeOutLength = .005;
  #endif
  std::string poissonType = "hypre";
  // output
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0;
  std::string path4serialization = "./";

  #ifdef CUP_ASYNC_DUMP
    MPI_Comm dump_comm;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
  #endif

  void resetAll();
  bool bDump();
  void registerDump();
  bool bOver() const;

  double minRho() const;
  double maxSpeed() const;
  double maxRelSpeed() const;
  void checkVariableDensity();

  inline double getH() const
  {
    return vel->getBlocksInfo().front().h_gridpoint; // yikes
  }

  void startProfiler(std::string name);
  void stopProfiler();
  void printResetProfiler();
  ~SimulationData();

  void dumpChi  (std::string name);
  void dumpPres (std::string name);
  void dumpPrhs (std::string name);
  void dumpTmp  (std::string name);
  void dumpVel  (std::string name);
  void dumpUobj (std::string name);
  void dumpForce(std::string name);
  void dumpTmpV (std::string name);
  void dumpAll  (std::string name);
};
