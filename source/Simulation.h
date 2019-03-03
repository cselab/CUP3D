//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Simulation_Fluid_h
#define CubismUP_3D_Simulation_Fluid_h

#include "SimulationData.h"
// Forward declarations.
namespace cubism { class ArgumentParser; }

class IF3D_ObstacleOperator;

class Simulation
{
  //#ifdef _USE_ZLIB_
  //  SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, ChiStreamer> waveletdumper_grid;
  //#endif

public:

  SimulationData sim;
  ArgumentParser * const parser_ptr = nullptr;

  void _init(bool restart = false);

  void _serialize(const std::string append = std::string());
  void _deserialize();

  void _argumentsSanityCheck();
  void setObstacleVector(IF3D_ObstacleVector *obstacle_vector_);
  void setupOperators();
  void setupGrid();
  void _ic();

 public:
  Simulation(const SimulationData &);
  Simulation(MPI_Comm mpicomm);
  Simulation(MPI_Comm mpicomm, cubism::ArgumentParser &parser);

  virtual ~Simulation() = default;

  virtual void run();

  // void addObstacle(IF3D_ObstacleOperator *obstacle);
  // void removeObstacle(IF3D_ObstacleOperator *obstacle);

  /* Get reference to the obstacle container. */
  const std::vector<IF3D_ObstacleOperator *> &getObstacleVector() const;

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
