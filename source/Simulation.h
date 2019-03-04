//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Written by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Simulation_h
#define CubismUP_3D_Simulation_h

#include "SimulationData.h"
// Forward declarations.
namespace cubism { class ArgumentParser; }
namespace cubismup3d { class SimulationWrapper; }

CubismUP_3D_NAMESPACE_BEGIN

class Obstacle;

class Simulation
{
  friend class cubismup3d::SimulationWrapper;
  //#ifdef _USE_ZLIB_
  //  SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, ChiStreamer> waveletdumper_grid;
  //#endif

public:

  SimulationData sim;
  cubism::ArgumentParser * const parser_ptr;

  void _init(bool restart = false);

  void _serialize(const std::string append = std::string());
  void _deserialize();

  void _argumentsSanityCheck();
  void setObstacleVector(ObstacleVector *obstacle_vector_);
  void setupOperators();
  void setupGrid();
  void _ic();

 public:
  Simulation(MPI_Comm mpicomm);
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
             double fadeOutLength,
             int saveFreq, double saveTime,
             const std::string &path4serialization,
             bool restart);

  virtual ~Simulation() {}

  virtual void run();

  // void addObstacle(IF3D_ObstacleOperator *obstacle);
  // void removeObstacle(IF3D_ObstacleOperator *obstacle);

  /* Get reference to the obstacle container. */
  const std::vector<Obstacle *> &getObstacleVector() const;

  /* Calculate maximum allowed time step, including CFL and ramp-up. */
  double calcMaxTimestep();

  /*
   * Perform one timestep of the simulation.
   *
   * Returns true if the simulation is finished.
   */
  bool timestep(double dt);
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Simulation_h
