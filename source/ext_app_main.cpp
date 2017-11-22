//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;
#include "Simulation.h"

int app_main(Communicator*const rlcom, MPI_Comm mpicom, int argc, char**argv)
{
  int rank;
  MPI_Comm_rank(mpicom, &rank);

  ArgumentParser parser(argc,argv);
  parser.set_strict_mode();

  int supported_threads;
  MPI_Query_thread(&supported_threads);
  if (supported_threads < MPI_THREAD_FUNNELED) {
    printf("ERROR: The MPI implementation does not have required thread support\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  printf("Rank %d is on host %s\n", rank, hostname);

  if (rank==0) {
   cout <<
    "=======================================================================\n";
   cout <<
    "Cubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
   cout <<
    "=======================================================================\n";
  }

        TaskLayer task(rlcom, mpicom, parser);
  Simulation sim(mpicom, &task, argc, argv);
  sim.init();
  sim.simulate();

  return 0;
}
