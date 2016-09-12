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

#include "Definitions.h"
#include "Simulation.h"

int main(int argc, const char **argv)
{
	//MPI_Init(&argc, &argv);
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if (rank==0) {
		cout << "====================================================================================================================\n";
		cout << "\t\tCubism UP 3D (velocity-pressure 3D incompressible Navier-Stokes solver)\n";
		cout << "====================================================================================================================\n";
	}
	
	ArgumentParser parser(argc,argv);
	parser.set_strict_mode();
	
	string simSetting = parser("-sim").asString();
	Simulation_Fluid * sim;
	
	sim = new Simulation(argc, argv);
	sim->init();
	sim->simulate();
	
	MPI_Finalize();
	
	return 0;
}
