//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_2D_Simulation_Fluid_h
#define CubismUP_2D_Simulation_Fluid_h

//#include "Definitions.h"
#include "GenericOperator.h"
#include "GenericCoordinator.h"

#include "IF3D_ObstacleVector.h"

#ifdef _USE_ZLIB_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <thread>
#include <vector>

#ifdef RL_LAYER
 #include "TaskLayer.h"
#endif

namespace cubismup3d { class SimulationWrapper; }

class Simulation
{
  friend class cubismup3d::SimulationWrapper;
 protected:
  ArgumentParser parser;
  Profiler profiler;
  const MPI_Comm app_comm;
  #ifdef RL_LAYER
    TaskLayer * task;
  #endif
  #ifdef _USE_ZLIB_
    SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, ChiStreamer> waveletdumper_grid;
  #endif

  // grid
  int rank=-1, nprocs=1;
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;

  // simulation status
  int step=0, nsteps=0;
  double time=0, endTime=0;

  // simulation settings
  Real uinf[3]={0,0,0};
  double nu=0, CFL=0, lambda=0, DLM=0;
  bool bRestart=false, verbose=false;
  bool b3Ddump=true, b2Ddump=false, bDump=false;
  bool rampup=true;

  // output
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0, saveClockPeriod=0, maxClockDuration=1e9;
  string path2file, path4serialization = "./";

  FluidGridMPI * grid = nullptr;

  #ifdef DUMPGRID
    MPI_Comm dump_comm;
    DumpGridMPI * dump = nullptr;
    std::thread * dumper = nullptr;
  #endif

  std::vector<BlockInfo> vInfo;
  //The protagonist
  IF3D_ObstacleVector* obstacle_vector;
  //The antagonist
  std::vector<GenericCoordinator*> pipeline;
  // vector of 2D slices (for dumping)
  std::vector<SliceType> m_slices;

  void _serialize(const string append = string());
  void _deserialize();

  void parseArguments();
  void setupObstacles();
  void setupOperators();
  void setupGrid();
  void _ic();

 public:
  Simulation(MPI_Comm mpicomm, int argc, char** argv):
  parser(argc,argv), app_comm(mpicomm) {   }

  virtual ~Simulation()
  {
    delete grid;
    delete obstacle_vector;
    while(!pipeline.empty()) {
      GenericCoordinator * g = pipeline.back();
      pipeline.pop_back();
      delete g;
    }
    #ifdef DUMPGRID
      if(dumper not_eq nullptr) {
        dumper->join();
        delete dumper;
      }
      delete dump;
      MPI_Comm_free(&dump_comm);
    #endif
  }

  virtual void init();
  virtual void simulate();


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
