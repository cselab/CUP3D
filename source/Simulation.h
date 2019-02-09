//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_2D_Simulation_Fluid_h
#define CubismUP_2D_Simulation_Fluid_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"

#ifdef _USE_ZLIB_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif

#include "Cubism/HDF5SliceDumperMPI.h"
#include "Cubism/Profiler.h"
//#include "Cubism/ZBinDumper_MPI.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <array>
#include <thread>
#include <vector>

#ifdef RL_LAYER
 #include "TaskLayer.h"
#endif

#ifdef CUP_ASYNC_DUMP
  using DumpBlock = BaseBlock<DumpElement>;
  typedef GridMPI<Grid<DumpBlock, aligned_allocator>> DumpGridMPI;
  typedef SliceTypesMPI::Slice<DumpGridMPI> SliceType;
#else
  typedef SliceTypesMPI::Slice<FluidGridMPI> SliceType;
#endif

// Forward declarations.
namespace cubism { class ArgumentParser; }
namespace cubismup3d { class SimulationWrapper; }

class Simulation
{
  friend class cubismup3d::SimulationWrapper;
 protected:
  Profiler profiler;
  const MPI_Comm app_comm;
  #ifdef RL_LAYER
    TaskLayer * task;
  #endif
  #ifdef _USE_ZLIB_
    SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, ChiStreamer> waveletdumper_grid;
  #endif

  int rank=-1, nprocs=-1;
public:
  // grid
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;

  // simulation status
  int step=0, nsteps=0;
  double time=0, endTime=0;

  // simulation settings
  Real uinf[3]={0,0,0};
  double nu=0, CFL=0, lambda=-1, DLM=1;
  bool verbose=false;
  bool computeDissipation=false;
  bool b3Ddump=true, b2Ddump=false, bDump=false;
  bool rampup=true;
#ifndef CUP_UNBOUNDED_FFT
  Real fadeOutLength = .005;
#endif

  // output
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0, saveClockPeriod=0, maxClockDuration=1e9;
  std::string path4serialization = "./";

  FluidGridMPI * grid = nullptr;

  #ifdef CUP_ASYNC_DUMP
    MPI_Comm dump_comm;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
  #endif

  std::vector<BlockInfo> vInfo;
  //The protagonist
  IF3D_ObstacleVector* obstacle_vector = nullptr;
  //The antagonist
  std::vector<GenericCoordinator*> pipeline;
  // vector of 2D slices (for dumping)
  std::vector<SliceType> m_slices;

  void _init(bool restart = false);

  void _serialize(const std::string append = std::string());
  void _deserialize();

  void _argumentsSanityCheck();
  void setObstacleVector(IF3D_ObstacleVector *obstacle_vector_);
  void setupOperators();
  void setupGrid();
  void _ic();

 public:
  Simulation(MPI_Comm mpicomm) : app_comm(mpicomm) {}
  Simulation(MPI_Comm mpicomm, cubism::ArgumentParser &parser);

  // For Python bindings. Order should be the same as defined in the class. The
  // default values are set in Python bindings.
  Simulation(std::array<int, 3> cells,
             std::array<int, 3> nproc,
             MPI_Comm mpicomm,
             int nsteps, double endTime,
             double nu, double CFL, double lambda, double DLM,
             std::array<double, 3> uinf,
             bool verbose,
             bool computeDissipation,
             bool b3Ddump, bool b2Ddump,
#ifndef CUP_UNBOUNDED_FFT
             double fadeOutLength,
#endif
             int saveFreq, double saveTime,
             const std::string &path4serialization,
             bool restart);

  virtual ~Simulation();

  virtual void run();

  // void addObstacle(IF3D_ObstacleOperator *obstacle);
  // void removeObstacle(IF3D_ObstacleOperator *obstacle);

  /* Get reference to the obstacle container. */
  const std::vector<IF3D_ObstacleOperator *> &getObstacleVector() const
  {
      return obstacle_vector->getObstacleVector();
  }

  /* Calculate maximum allowed time step, including CFL and ramp-up. */
  double calcMaxTimestep();

  /*
   * Perform one timestep of the simulation.
   *
   * Returns true if the simulation is finished.
   */
  bool timestep(double dt);
};

#endif
