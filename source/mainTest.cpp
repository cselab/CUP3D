//
//  mainTest.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
using namespace std;

#include "Definitions.h"

#include "TestDiffusion.h"
#include "TestAdvection.h"
#include "TestPressure.h"
#include "TestVarCoeffPoisson.h"
#include "TestGravity.h"
#include "TestPenalization.h"
#include "TestTranslation.h"
#include "TestRotation.h"
#include "TestTravelingWave.h"
#include "TestShearLayer.h"
#include "TestPoiseuille.h"
#include "TestAddedMass.h"
#include "TestBoundaryConditions.h"
#include "TestGeometry.h"
#include "TestMPI.h"

void spatialConvergence(int argc, const char **argv, const int solver, const int ic, const string test, const int minBPD, const int maxBPD, const double dt)
{
	if (test=="advection")
	{
		if (ic==0)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test (Space) - Linear Field\n";
			cout << "========================================================================================\n";
		}
		else if (ic==1)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test (Space) - Vortex Field\n";
			cout << "========================================================================================\n";
		}
		else
			throw std::invalid_argument("This test setting does not exist!");
		
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestAdvection * advection = new TestAdvection(argc, argv, ic, bpd, dt, 1);
			advection->run();
			advection->check();
			delete advection;
		}
	}
	else if (test=="diffusion")
	{
		cout << "========================================================================================\n";
		cout << "\t\tDiffusion Test (Space)\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestDiffusion * diffusion = new TestDiffusion(argc, argv, bpd, dt, 1, 8);
			diffusion->run();
			diffusion->check();
			delete diffusion;
		}
	}
	else if (test=="pressure")
	{
		cout << "========================================================================================\n";
		if (solver==0)
			cout << "\t\tPressure Test (Space) - Stencil - ";
		else if (solver==1)
			cout << "\t\tPressure Test (Space) - Spectral - ";
#ifdef _MULTIGRID_
		else if (solver==2)
			cout << "\t\tPressure Test (Space) - Multigrid (Constant Coefficents) - ";
		else if (solver==3)
			cout << "\t\tPressure Test (Space) - Multigrid (Variable Coefficents) - ";
#endif // _MULTIGRID_
		else
			throw std::invalid_argument("This test setting does not exist!");
		
		if (ic==0)
			cout << "Poisson equation\n";
		else if (ic==1)
			cout << "Velocity field (Projection)\n";
		else if (ic==2)
			cout << "Mixed Boundary Conditions\n";
		else
			throw std::invalid_argument("This test setting does not exist!");
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestPressure * pressure = new TestPressure(argc, argv, solver, ic, bpd, dt);
			pressure->run();
			pressure->check();
			delete pressure;
		}
	}
	else if (test=="poisson")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPoisson Test (Space) - Multigrid - ";
		
		if (ic==0)
			cout << "Constant Coefficients\n";
		else if (ic==1)
			cout << "Variable Coefficients\n";
		else
			throw std::invalid_argument("This test setting does not exist!");
		cout << "========================================================================================\n";
		
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestVarCoeffPoisson * poisson = new TestVarCoeffPoisson(argc, argv, ic, bpd);
			poisson->run();
			poisson->check();
			delete poisson;
		}
	}
	else if (test=="translation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tTranslation Test (Space) - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tTranslation Test (Space) - Velocity from Flow\n";
		else
			throw std::invalid_argument("This test setting does not exist!");
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestTranslation * translation = new TestTranslation(argc, argv, ic, bpd, dt);
			translation->run();
			translation->check();
			delete translation;
		}
	}
	else if (test=="rotation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tRotation Test (Space) - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tRotation Test (Space) - Velocity from Flow\n";
		else
			throw std::invalid_argument("This test setting does not exist!");
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestRotation * rotation = new TestRotation(argc, argv, ic, bpd, dt);
			rotation->run();
			rotation->check();
			delete rotation;
		}
	}
	else if (test=="travelingwave")
	{
		cout << "========================================================================================\n";
		cout << "\t\tTraveling Wave Test (Space)\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestTravelingWave * wave = new TestTravelingWave(argc, argv, bpd);
			wave->run();
			wave->check();
			delete wave;
		}
	}
	else if (test=="addedmass")
	{
		cout << "========================================================================================\n";
		cout << "\t\tAdded Mass Test (Space)\n";
		cout << "========================================================================================\n";
		for (int bpd=minBPD; bpd<=maxBPD; bpd*=2)
		{
			TestAddedMass * am = new TestAddedMass(argc, argv, bpd);
			am->run();
			am->check();
			delete am;
		}
	}
	else
		throw std::invalid_argument("This test setting does not exist!");
}

void temporalConvergence(int argc, const char **argv, const int solver, const int ic, const string test, const double minDT, const double maxDT, const int bpd)
{
	if (test=="advection")
	{
		if (ic==0)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test (Time) - Linear Field\n";
			cout << "========================================================================================\n";
		}
		else if (ic==1)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test (Time) - Vortex Field\n";
			cout << "========================================================================================\n";
		}
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestAdvection * advection = new TestAdvection(argc, argv, ic, bpd, dt, 1000);
			advection->run();
			advection->check();
			delete advection;
		}
	}
	else if (test=="diffusion")
	{
		cout << "========================================================================================\n";
		cout << "\t\tDiffusion Test (Time)\n";
		cout << "========================================================================================\n";
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestDiffusion * diffusion = new TestDiffusion(argc, argv, bpd, dt, 10000, 4);
			diffusion->run();
			diffusion->check();
			delete diffusion;
		}
	}
	else if (test=="gravity")
	{
		cout << "========================================================================================\n";
		cout << "\t\tGravity Test (Time)\n";
		cout << "========================================================================================\n";
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestGravity * gravity = new TestGravity(argc, argv, bpd, dt);
			gravity->run();
			gravity->check();
			delete gravity;
		}
	}
	else if (test=="penalization")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPenalization Test (Time)\n";
		cout << "========================================================================================\n";
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestPenalization * penalization = new TestPenalization(argc, argv, bpd, dt);
			penalization->run();
			delete penalization;
		}
	}
	else if (test=="translation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tTranslation Test (Time) - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tTranslation Test (Time) - Velocity from Flow\n";
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		cout << "========================================================================================\n";
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestTranslation * translation = new TestTranslation(argc, argv, ic, bpd, dt);
			translation->run();
			translation->check();
			delete translation;
		}
	}
	else if (test=="rotation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tRotation Test (Time) - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tRotation Test (Time) - Velocity from Flow\n";
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		cout << "========================================================================================\n";
		for (double dt=minDT; dt<=maxDT; dt*=2)
		{
			TestRotation * rotation = new TestRotation(argc, argv, ic, bpd, dt);
			rotation->run();
			rotation->check();
			delete rotation;
		}
	}
	else
		throw std::invalid_argument("This test setting does not exist!");
}

void baseTest(int argc, const char **argv, const int solver, const int ic, const string test, const int bpd, const double dt, const double tEnd)
{
	if (test=="advection")
	{
		if (ic==0)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test - Linear Field\n";
			cout << "========================================================================================\n";
		}
		else if (ic==1)
		{
			cout << "========================================================================================\n";
			cout << "\t\tAdvection Test - Vortex Field\n";
			cout << "========================================================================================\n";
		}
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		
		TestAdvection * advection = new TestAdvection(argc, argv, ic, bpd, dt, tEnd);
		advection->run();
		advection->check();
		delete advection;
	}
	else if (test=="diffusion")
	{
		cout << "========================================================================================\n";
		cout << "\t\tDiffusion Test\n";
		cout << "========================================================================================\n";
		TestDiffusion * diffusion = new TestDiffusion(argc, argv, bpd, dt, tEnd, 8);
		diffusion->run();
		diffusion->check();
		delete diffusion;
	}
	else if (test=="pressure")
	{
		cout << "========================================================================================\n";
		if (solver==0)
			cout << "\t\tPressure Test - Stencil - ";
		else if (solver==1)
			cout << "\t\tPressure Test - Spectral - ";
#ifdef _MULTIGRID_
		else if (solver==2)
			cout << "\t\tPressure Test - Multigrid (Constant Coefficents) - ";
		else if (solver==3)
			cout << "\t\tPressure Test - Multigrid (Variable Coefficents) - ";
#endif // _MULTIGRID_
		else
			throw std::invalid_argument("Chosen solver does not exist!");
		
		if (ic==0)
			cout << "Poisson equation\n";
		else if (ic==1)
			cout << "Velocity field (Projection)\n";
		else if (ic==2)
			cout << "Mixed Boundary Conditions\n";
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		cout << "========================================================================================\n";
		TestPressure * pressure = new TestPressure(argc, argv, solver, ic, bpd, dt);
		pressure->run();
		pressure->check();
		delete pressure;
	}
	else if (test=="gravity")
	{
		cout << "========================================================================================\n";
		cout << "\t\tGravity Test\n";
		cout << "========================================================================================\n";
		TestGravity * gravity = new TestGravity(argc, argv, bpd, dt);
		gravity->run();
		gravity->check();
		delete gravity;
	}
	else if (test=="penalization")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPenalization Test\n";
		cout << "========================================================================================\n";
		TestPenalization * penalization = new TestPenalization(argc, argv, bpd, dt);
		penalization->run();
		delete penalization;
	}
	else if (test=="translation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tTranslation Test - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tTranslation Test - Velocity from Flow\n";
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		cout << "========================================================================================\n";
		TestTranslation * translation = new TestTranslation(argc, argv, ic, bpd, dt);
		translation->run();
		translation->check();
		delete translation;
	}
	else if (test=="rotation")
	{
		cout << "========================================================================================\n";
		if (ic==0)
			cout << "\t\tRotation Test - Constant Velocity\n";
		else if (ic==1)
			cout << "\t\tRotation Test - Velocity from Flow\n";
		else
			throw std::invalid_argument("Chosen IC does not exist!");
		cout << "========================================================================================\n";
		TestRotation * rotation = new TestRotation(argc, argv, ic, bpd, dt);
		rotation->run();
		rotation->check();
		delete rotation;
	}
	else if (test=="travelingwave")
	{
		cout << "========================================================================================\n";
		cout << "\t\tTraveling Wave Test\n";
		cout << "========================================================================================\n";
		TestTravelingWave * wave = new TestTravelingWave(argc, argv, bpd);
		wave->run();
		wave->check();
		delete wave;
	}
	else if (test=="shearlayer")
	{
		cout << "========================================================================================\n";
		cout << "\t\tShear Layer Test\n";
		cout << "========================================================================================\n";
		TestShearLayer * shearlayer = new TestShearLayer(argc, argv, bpd);
		shearlayer->run();
		shearlayer->check();
		delete shearlayer;
	}
	else if (test=="poiseuille")
	{
		cout << "========================================================================================\n";
		cout << "\t\tPoiseuille Test\n";
		cout << "========================================================================================\n";
		TestPoiseuille * poiseuille = new TestPoiseuille(argc, argv, bpd);
		poiseuille->run();
		poiseuille->check();
		delete poiseuille;
    }
    else if (test=="addedmass")
    {
        cout << "========================================================================================\n";
        cout << "\t\tAdded Mass Test\n";
        cout << "========================================================================================\n";
        TestAddedMass * am = new TestAddedMass(argc, argv, bpd);
        am->run();
        am->check();
        delete am;
	}
	else if (test=="bc")
	{
		cout << "========================================================================================\n";
		cout << "\t\tBoundary Conditions Test\n";
		cout << "========================================================================================\n";
		TestBoundaryConditions * am = new TestBoundaryConditions(argc, argv);
		am->run();
		am->check();
		delete am;
	}
	else if (test=="geometry")
	{
		cout << "========================================================================================\n";
		cout << "\t\tGeometry Test\n";
		cout << "========================================================================================\n";
		TestGeometry * geometry = new TestGeometry(argc, argv, bpd);
		geometry->run();
		delete geometry;
	}
	else if (test=="mpi")
	{
		cout << "========================================================================================\n";
		cout << "\t\tMPI Test\n";
		cout << "========================================================================================\n";
		TestMPI * mpi = new TestMPI(argc, argv, bpd);
		mpi->run();
		mpi->check();
		delete mpi;
	}
	else
		throw std::invalid_argument("This test setting does not exist!");
}

int main(int argc, const char **argv)
{
	MPI_Init(&argc, &argv);
	
	ArgumentParser parser(argc,argv);
	int solver = parser("-solver").asInt(0);
	int ic = parser("-ic").asInt(0);
	const double tEnd = parser("-tEnd").asDouble(0);
	
	const double minDT = parser("-minDT").asDouble(0);
	const double maxDT = parser("-maxDT").asDouble(0);
	parser.set_strict_mode();
	
	string test = parser("-test").asString();
	const int minBPD = parser("-minBPD").asInt();
	const int maxBPD = parser("-maxBPD").asInt();
	
	if (maxBPD>minBPD)
		spatialConvergence(argc, argv, solver, ic, test, minBPD, maxBPD, minDT);
	else if (maxDT>minDT)
		temporalConvergence(argc, argv, solver, ic, test, minDT, maxDT, minBPD);
	else
		baseTest(argc, argv, solver, ic, test, minBPD, minDT, tEnd);
	
	MPI_Finalize();
	
	return 0;
}