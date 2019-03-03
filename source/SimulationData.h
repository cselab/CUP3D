//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#pragma once
#include "Definitions.h"
#include "Cubism/Profiler.h"
#include "utils/NonUniformScheme.h"
#include "Cubism/HDF5SliceDumperMPI.h"
//#include "Cubism/ZBinDumper_MPI.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <array>
#include <thread>
#include <vector>

class Operator;
class IF3D_ObstacleVector;
class PoissonSolver;

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
  NonUniformScheme<FluidBlock>* nonuniform = nullptr;
  const std::vector<BlockInfo>& vInfo() const {
    return grid->getBlocksInfo();
  }

  // vector of 2D slices (for dumping)
  std::vector<SliceType> m_slices;

  //The protagonist
  IF3D_ObstacleVector * obstacle_vector = nullptr;
  //The antagonist
  std::vector<Operator*> pipeline;
  PoissonSolver * pressureSolver = nullptr;
  // simulation status
  // nsteps==0 means that this stopping criteria is not active
  int step=0, nsteps=0;
  // endTime==0  means that this stopping criteria is not active
  double time=0, endTime=0;
  double dt = 0;

  // mpi
  const MPI_Comm app_comm;
  int rank=-1, nprocs=-1;
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;

  // grid
  int local_bpdx=-1, local_bpdy=-1, local_bpdz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;
  Real maxextent = 1;
  Real extent[3] = {1, 1, 1};
  bool bUseStretchedGrid = false;
  Real hmin=0;

  // flow variables
  Real uinf[3] = {0, 0, 0};
  double nu=0, CFL=0, lambda=-1, DLM=1;

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
  std::string useSolver = "";
  // flags assume value 0 for dirichlet/unbounded, 1 for periodic, 2 for wall
  BCflag BCx_flag = dirichlet, BCy_flag = dirichlet, BCz_flag = dirichlet;

  bool bUseUnboundedBC = false;
  bool bUseFourierBC = false;

  #ifdef CUP_ASYNC_DUMP
    MPI_Comm dump_comm = MPI_COMM_NULL;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
    NonUniformScheme<DumpBlock>* nonuniform_dump = nullptr;
  #endif

  void startProfiler(std::string name) const;
  void stopProfiler() const;
  void printResetProfiler();
  void _preprocessArguments();
  ~SimulationData();
  SimulationData(MPI_Comm mpicomm, ArgumentParser &parser);
  SimulationData(MPI_Comm mpicomm);
  void setCells(int nx, int ny, int nz);
};
