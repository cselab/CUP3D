//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#pragma once

#include "operators/GenericOperator.h"
#include "Cubism/Profiler.h"
#include "Cubism/HDF5SliceDumperMPI.h"
//#include "Cubism/ZBinDumper_MPI.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <array>
#include <thread>
#include <vector>

class GenericCoordinator;
class IF3D_ObstacleVector;

#ifdef _USE_ZLIB_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif

#ifdef CUP_ASYNC_DUMP
  using DumpBlock = BaseBlock<DumpElement>;
  typedef GridMPI<Grid<DumpBlock, aligned_allocator>> DumpGridMPI;
  typedef SliceTypesMPI::Slice<DumpGridMPI> SliceType;
#else
  typedef SliceTypesMPI::Slice<FluidGridMPI> SliceType;
#endif

struct SimulationData
{
  #ifndef SMARTIES_APP
    Profiler * profiler = new Profiler();
  #endif

  FluidGridMPI * grid = nullptr;
  const std::vector<BlockInfo>& vInfo() const {
    return grid->getBlocksInfo();
  }

  // vector of 2D slices (for dumping)
  std::vector<SliceType> m_slices;

  //The protagonist
  IF3D_ObstacleVector * obstacle_vector = nullptr;
  //The antagonist
  std::vector<GenericCoordinator*> pipeline;

  // simulation status
  // nsteps==0 means that this stopping criteria is not active
  int step=0, nsteps=0;
  // endTime==0  means that this stopping criteria is not active
  double time=0, endTime=0;
  double dt = 0;
  // mpi
  const MPI_Comm app_comm;
  int rank=-1, nprocs=-1;
  // grid
  int bpdx=-1, bpdy=-1, bpdz=-1;
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;
  // flow variables
  Real uinf[3] = {0, 0, 0};
  double nu=0, CFL=0, lambda=-1, DLM=1;
  const Real maxextent = 1;//grid->maxextent; TODO
  Real extent[3] = {1, 1, 1};
  // simulation settings
  bool computeDissipation=false;
  bool b3Ddump=true, b2Ddump=false, bDump=false;
  int rampup = 100;
  bool verbose=false;
  bool muteAll = false;
  Real fadeOutLengthU[3] = {0, 0, 0};
  Real fadeOutLengthPRHS[3] = {0, 0, 0};
  Real uMax_forced = 0;

  // output
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0;
  std::string path4serialization = "./";
  std::string initCond = "zero";
  // flags assume value 0 for dirichlet/unbounded, 1 for periodic, 2 for wall
  BCflag BCx_flag = dirichlet, BCy_flag = dirichlet, BCz_flag = dirichlet;

  bool bUseUnboundedBC = false;
  bool bUseFourierBC = false;

  #ifdef CUP_ASYNC_DUMP
    MPI_Comm dump_comm;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
  #endif

  void startProfiler(std::string name);
  void stopProfiler();
  void printResetProfiler();
  void _argumentsSanityCheck();
  ~SimulationData();
  SimulationData(MPI_Comm mpicomm, ArgumentParser &parser);
  SimulationData(MPI_Comm mpicomm);
};
