//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include <Cubism/ArgumentParser.h>
#include "SimulationData.h"

CubismUP_3D_NAMESPACE_BEGIN

class Simulation
{
protected:
  cubism::ArgumentParser parser;
public:
  SimulationData sim;

  void initialGridRefinement();
  void serialize(const std::string append = std::string());
  void deserialize();
  void setupOperators();
  void setupGrid();
  void _ic();

  //Simulation(MPI_Comm mpicomm, cubism::ArgumentParser &parser);
  Simulation(int argc, char ** argv, MPI_Comm comm);

  void init();

  void simulate();

  /// Manually trigger mesh adaptation.
  void adaptMesh();

  /* Get reference to the obstacle container. */
  const std::vector<std::shared_ptr<Obstacle>> &getShapes() const;

  /* Calculate maximum allowed time step, including CFL and ramp-up. */
  Real calcMaxTimestep();

  /*
   * Perform one timestep of the simulation.
   *
   * Returns true if the simulation is finished.
   */
  bool advance(Real dt);

  /// Compute vorticity and store to tmpU, tmpV and tmpW.
  void computeVorticity();

  /// Insert the operator at the end of the pipeline.
  void insertOperator(std::shared_ptr<Operator> op);
};

/** Create a Simulation object from a vector of command-line arguments.

    The argv vector should NOT contain the argv[0] argument, it is filled with
    a dummy value instead.
*/
std::shared_ptr<Simulation> createSimulation(
    MPI_Comm comm,
    const std::vector<std::string> &argv);

CubismUP_3D_NAMESPACE_END
