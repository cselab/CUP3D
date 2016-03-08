//
//  TestAddedMass.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 6/18/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestAddedMass.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"
#include "CoordinatorBodyVelocities.h"

void TestAddedMass::_ic()
{
	bSplit = parser("-split").asBool(false);
	nu = parser("-nu").asDouble(1e-2);
	lambda = parser("-lambda").asDouble(1e5);
	path2file = parser("-file").asString("../data/AddedMass");
	
	rhoS = parser("-rhoS").asDouble(1);
	minRho = min((Real)1.,(Real)rhoS);
	
	Real centerOfMass[3] = {.5,.5,.5};
	
    Real radius = parser("-radius").asDouble(0.1);
	shape = new Sphere(centerOfMass, radius, rhoS, 2, 2);
	
	if (rank==0)
	{
		// setup initial conditions
		CoordinatorIC coordIC(shape,0,grid);
		coordIC(0);
	}
}

TestAddedMass::TestAddedMass(const int argc, const char ** argv, const int bpd) : Test(argc, argv, bpd), parser(argc,argv), gravity{0,-9.81,0}, uBody{0,0,0}, bSplit(false), step(0), rank(0), nprocs(1), maxU(0)
{
	// output settings
	path2file = parser("-file").asString("../data/TestAddedMass");
	
	_ic();
	
	pipeline.clear();
	
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0], &uBody[1], &uBody[2], grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0], &uBody[1], &uBody[2], grid,1));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, gravity, &step, bSplit, grid, rank, nprocs));
	pipeline.push_back(new CoordinatorBodyVelocities(&uBody[0], &uBody[1], &uBody[2], &lambda, shape, &maxU, grid));
	pipeline.push_back(new CoordinatorComputeShape(shape, grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &uBody[2], shape, &lambda, grid));
}

TestAddedMass::~TestAddedMass()
{
	while(!pipeline.empty())
	{
		GenericCoordinator * g = pipeline.back();
		pipeline.pop_back();
		delete g;
	}
}

void TestAddedMass::run()
{
	double vOld = 0;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	nsteps = parser("-nsteps").asInt(3);
	double dt = 1e-07;
	
	for (step=0; step<nsteps; step++)
	{
		for (int c=0; c<pipeline.size(); c++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			
			if (rank == 0 || pipeline[c]->getName()=="Pressure")
				(*pipeline[c])(dt);
		}
		
		if (rank==0 && step==nsteps-1)
		{
#ifdef _USE_HDF_
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream sstmp;
			sstmp << path2file << "-" << std::setfill('0') << std::setw(6) << step;
			cout << sstmp.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, sstmp.str());
#endif
			
			// this still needs to be corrected to the frame of reference!
			double accM = (uBody[1]-vOld)/dt;
			vOld = uBody[1];
			double accT = (shape->getRhoS()-1)/(shape->getRhoS()+.5) * gravity[1];
			double accN = (shape->getRhoS()-1)/(shape->getRhoS()   ) * gravity[1];
			cout << bpd << "\t" << rhoS << "\t" << step << "\t" << accM << "\t" << accT << "\t" << accN << endl;
			stringstream ss;
			ss << path2file << "_addedmass.dat";
			ofstream myfile(ss.str(), fstream::app);
			myfile << bpd << "\t" << rhoS << "\t" << step << "\t" << accM << "\t" << accT << "\t" << accN << endl;
		}
	}
}

void TestAddedMass::check()
{
}