//
//  main.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;
#define __SMARTIES_
#include "Simulation.h"

int main(int argc, char **argv)
{
	Communicator * communicator = nullptr;
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if (provided < MPI_THREAD_FUNNELED) {
		printf("ERROR: The MPI implementation does not have required thread support\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	ArgumentParser parser(argc,argv);
	parser.set_strict_mode();

	if (rank==0) {
		cout << "====================================================================================================================\n";
		cout << "\t\tCubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
		cout << "====================================================================================================================\n";
#ifdef __SMARTIES_
		parser.unset_strict_mode();
		const int _sockID = parser("-sock").asInt(-1);
		const int nActions = parser("-nActions").asInt(0);
		const int nStates = (nActions==1) ? 20+200 : 25+200;
		if (_sockID>=0 && nActions>0) {
			printf("Communicating over sock %d\n", _sockID);
			communicator = new Communicator(_sockID,nStates,nActions);
		}
#endif
	}

	Simulation * sim = new Simulation(argc, argv, communicator);
	sim->init();
	sim->simulate();

	MPI_Finalize();

	return 0;
}
