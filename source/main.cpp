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
#include "Save_splicer.h"

int main(int argc, char **argv)
{
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) {
    printf("ERROR: MPI implementation does not have required thread support\n");
    fflush(0);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  ArgumentParser parser(argc,argv);
  parser.set_strict_mode();
  if (rank==0) {
   cout <<
    "=======================================================================\n";
   cout <<
    "Cubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
   cout <<
    "=======================================================================\n";
  }
  parser.unset_strict_mode();

  Simulation * sim = new Simulation(MPI_COMM_WORLD, argc, argv);
  sim->init();
  sim->simulate();

  MPI_Finalize();
  delete sim;
  return 0;
}
