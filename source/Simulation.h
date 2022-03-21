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

#include <memory>

// Forward declarations.
namespace cubism { class ArgumentParser; }

CubismUP_3D_NAMESPACE_BEGIN

class Obstacle;

class Simulation
{
public:

  SimulationData sim;

  void initialGridRefinement();
  void _serialize(const std::string append = std::string());
  void _deserialize();

  void _argumentsSanityCheck();
  void setupOperators();
  void setupGrid();
  void _ic();
  void _icFromH5(std::string h5File);

  Simulation(MPI_Comm mpicomm, cubism::ArgumentParser &parser);
  Simulation(const Simulation &) = delete;
  Simulation(Simulation &&) = delete;
  Simulation &operator=(const Simulation &) = delete;
  Simulation &operator=(Simulation &&) = delete;

  virtual ~Simulation() = default;

  virtual void run();

  /// Manually trigger mesh adaptation.
  void adaptMesh();

  /* Get reference to the obstacle container. */
  const std::vector<std::shared_ptr<Obstacle>> &getObstacleVector() const;

  /* Calculate maximum allowed time step, including CFL and ramp-up. */
  Real calcMaxTimestep();

  /*
   * Perform one timestep of the simulation.
   *
   * Returns true if the simulation is finished.
   */
  bool timestep(Real dt);

  /// Insert the operator at the end of the pipeline.
  void insertOperator(std::shared_ptr<Operator> op);

  void touch();
};

/** Create a Simulation object from a vector of command-line arguments.

    The argv vector should NOT contain the argv[0] argument, it is filled with
    a dummy value instead.
*/
std::shared_ptr<Simulation> createSimulation(
    MPI_Comm comm,
    const std::vector<std::string> &argv);

CubismUP_3D_NAMESPACE_END

#endif // CubismUP_3D_Simulation_h
