//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_SimulationData_h
#define CubismUP_3D_SimulationData_h

#include "Definitions.h"
#ifdef _USE_ZLIB_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif
//#include "Cubism/ZBinDumper_MPI.h"

#include <array>
#ifdef CUP_ASYNC_DUMP
#include <thread>
#endif
#include <vector>

namespace cubism {
  class Profiler;
  class ArgumentParser;
  namespace SliceTypesMPI { template<typename grid_t> class Slice; }
}

CubismUP_3D_NAMESPACE_BEGIN

class Operator;
class ObstacleVector;
class PoissonSolver;

#ifdef CUP_ASYNC_DUMP
 using DumpBlock  = BaseBlock<DumpElement>;
 using DumpGridMPI= cubism::GridMPI<cubism::Grid<DumpBlock, aligned_allocator>>;
 using SliceType  = cubism::SliceTypesMPI::Slice<DumpGridMPI>;
#else
 using SliceType  = cubism::SliceTypesMPI::Slice<FluidGridMPI>;
#endif

struct SimulationData
{
  cubism::Profiler * profiler = nullptr;

  FluidGridMPI * grid = nullptr;
  void * nonuniform = nullptr;
  const inline std::vector<cubism::BlockInfo>& vInfo() const {
    return grid->getBlocksInfo();
  }
  Real maxH() const { return hmax; }
  Real uniformH() const
  {
    if(std::fabs(hmin-hmax) > 1e-15) {
      printf("ABORT: SimulationData::uniformH requires uniform grids.\n");
      fflush(0); MPI_Abort(grid->getCartComm(), 1);
    }
    return hmin;
  }

  // vector of 2D slices (for dumping)
  std::vector<SliceType> m_slices;

  //The protagonist
  ObstacleVector * obstacle_vector = nullptr;
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
  MPI_Comm app_comm;
  int rank=-1, nprocs=-1;
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;

  // grid
  int local_bpdx=-1, local_bpdy=-1, local_bpdz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;
  Real maxextent = 1;
  std::array<Real, 3> extent = {{1, 0, 0}};  // Uniform grid by default.
  bool bUseStretchedGrid = false;
  bool bIterativePenalization = false;
  Real hmin=0, hmax=0;

  // flow variables
  std::array<Real, 3> uinf = {{0, 0, 0}};
  double nu=0, CFL=0, lambda=-1, DLM=1;

  // simulation settings
  int freqDiagnostics = 0;
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
    void * dump_nonuniform = nullptr;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
  #endif

  void startProfiler(std::string name) const;
  void stopProfiler() const;
  void printResetProfiler();
  void _preprocessArguments();
  ~SimulationData();
  SimulationData() = delete;
  SimulationData(const SimulationData &);
  SimulationData(SimulationData &&);
  SimulationData(MPI_Comm mpicomm, cubism::ArgumentParser &parser);
  SimulationData(MPI_Comm mpicomm);
  void setCells(int nx, int ny, int nz);
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SimulationData_h
