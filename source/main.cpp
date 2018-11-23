//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include <iostream>

#include "Cubism/ArgumentParser.h"
#include "IF3D_ObstacleFactory.h"
#include "Simulation.h"

int main(int argc, char **argv)
{
  int provided;
  #ifdef DUMPGRID
    const auto SECURITY = MPI_THREAD_MULTIPLE;
  #else
    const auto SECURITY = MPI_THREAD_FUNNELED;
  #endif
  MPI_Init_thread(&argc, &argv, SECURITY, &provided);
  if (provided < SECURITY ) {
    printf("ERROR: MPI implementation does not have required thread support\n");
    fflush(0);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank==0) {
    std::cout << "=======================================================================\n";
    std::cout << "Cubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
    std::cout << "=======================================================================\n";
#ifdef NDEBUG
    std::cout << "Running in RELEASE mode!\n";
#else
    std::cout << "Running in DEBUG mode!\n";
#endif
#ifdef _UNBOUNDED_FFT_
    std::cout << "Using freespace unbounded FFT...\n";
#else
    std::cout << "Using smooth truncation at domain boundaries...\n";
#endif
  }

  ArgumentParser parser(argc, argv);
  Simulation *sim = new Simulation(MPI_COMM_WORLD, parser);
  sim->run();
  delete sim;

  MPI_Finalize();
  return 0;
}
