//
//  Simulation_Fluid.h
//  CubismUP_2D
//
//  Base class for fluid simulations from which any fluid simulation case should inherit
//  Contains the base structure and interface that any fluid simulation class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_Simulation_Fluid_h
#define CubismUP_2D_Simulation_Fluid_h

//#include "Definitions.h"
#include "GenericOperator.h"
#include "GenericCoordinator.h"

#include "CoordinatorIC.h"
//#include "CoordinatorVorticity.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorFadeOut.h"
#include "IF3D_ObstacleVector.h"
#include "IF3D_ObstacleFactory.h"
//#include "IF3D_StefanFishOperator.h"
//#include "IF3D_CarlingFishOperator.h"
//#include "IF3D_SphereObstacleOperator.h"
//#include "IF3D_ForcedSphereObstacleOperator.h"
#if _USE_ZLIB_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <thread>

#ifdef RL_LAYER
 #include "TaskLayer.h"
#endif

class Simulation
{
 protected:
  ArgumentParser parser;
  Profiler profiler;
  const MPI_Comm app_comm;
  #ifdef RL_LAYER
    TaskLayer * task;
  #endif
  #if _USE_ZLIB_
    SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, ChiStreamer> waveletdumper_grid;
  #endif

  // grid
  int rank=-1, nprocs=1;
  int nprocsx=-1, nprocsy=-1, nprocsz=-1;
  int bpdx=-1, bpdy=-1, bpdz=-1;

  // simulation status
  int step=0, nsteps=0;
  double dt=0, time=0, endTime=0, dtCFL=0, dtFourier=0;

  // simulation settings
  Real uinf[3]={0,0,0};
  double re=0, nu=0, length=0, CFL=0, lambda=0, DLM=0;
  bool bRestart=false, verbose=false;
  bool b3Ddump=true, b2Ddump=false, bDump=false;

  // output
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0, saveClockPeriod=0, maxClockDuration=1e9;
  string path2file, path4serialization = "./";

  FluidGridMPI * grid = nullptr;
  DumpGridMPI * dump = nullptr;
  std::thread * dumper = nullptr;

  vector<BlockInfo> vInfo;
  //The protagonist
  IF3D_ObstacleVector* obstacle_vector;
  //The antagonist
  vector<GenericCoordinator*> pipeline;

  void _serialize(const string append = string());
  void _deserialize();

  void parseArguments();
  void setupObstacles();
  void setupOperators();
  void _selectDT();
  void setupGrid();
  void _ic();

 public:
  Simulation(MPI_Comm mpicomm, int argc, char** argv):
  parser(argc,argv), app_comm(mpicomm) {   }

  virtual ~Simulation()
  {
    delete grid;
    delete dump;
    while(!pipeline.empty()) {
      GenericCoordinator * g = pipeline.back();
      pipeline.pop_back();
      delete g;
    }
  }

  virtual void init();
  virtual void simulate();
};

#endif
