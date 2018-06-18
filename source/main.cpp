//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;
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

  ArgumentParser parser(argc,argv);
  parser.set_strict_mode();
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
#endif
  }
  parser.unset_strict_mode();

  Simulation * sim = new Simulation(MPI_COMM_WORLD, argc, argv);
  sim->init();
  sim->simulate();
  delete sim;

  MPI_Finalize();
  return 0;
}
