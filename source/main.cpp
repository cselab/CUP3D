//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Simulation.h"

#include <Cubism/ArgumentParser.h>

#include <iostream>
int main(int argc, char **argv)
{
  int provided;
  const auto SECURITY = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, SECURITY, &provided);
  if (provided < SECURITY ) {
    printf("ERROR: MPI implementation does not have required thread support\n");
    fflush(0); MPI_Abort(MPI_COMM_WORLD, 1);
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
  }
 
  MPI_Barrier(MPI_COMM_WORLD);
  double t1 = MPI_Wtime();
  cubism::ArgumentParser parser(argc, argv);
  cubismup3d::Simulation *sim = new cubismup3d::Simulation(MPI_COMM_WORLD, parser);
  sim->run();
  delete sim;
  MPI_Barrier(MPI_COMM_WORLD);
  double t2 = MPI_Wtime();
  if (rank == 0) std::cout << "Total time = " << t2 - t1 << std::endl;

  MPI_Finalize();
  return 0;
}
