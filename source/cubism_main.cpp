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
#include "Simulation.h"
#include "Save_splicer.h"

int cubism_main (const MPI_Comm app_comm, int argc, char **argv)
{
	Communicator * communicator = nullptr;

	int rank;
	MPI_Comm_rank(app_comm, &rank);

	ArgumentParser parser(argc,argv);
	parser.set_strict_mode();

	if (rank==0) {
		cout << "====================================================================================================================\n";
		cout << "\t\tCubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
		cout << "====================================================================================================================\n";
	}

	#ifdef __SMARTIES_
	parser.unset_strict_mode();
	const int _sockID = parser("-sock").asInt(-1);
	const int nActions = parser("-nActions").asInt(0);
	const int nStates = (nActions==1) ? 20+200 : 25+200;
	if (_sockID>=0 && nActions>0) {
		if(!rank)
			printf("Communicating over sock %d\n", _sockID);
		communicator = new Communicator(_sockID, nStates, nActions, app_comm);
	}
	#endif

	const bool io_job = parser("-saveSplicer").asBool(false);
	if (io_job) {
		Save_splicer * sim = new Save_splicer(argc, argv);
		sim->run();
		MPI_Finalize();
		return 0;
	}

	Simulation * sim = new Simulation(argc, argv, app_comm, communicator);
	sim->init();
	sim->simulate();

	MPI_Finalize();
	if(communicator not_eq nullptr)
			delete communicator;
	delete sim;
	return 0;
}
