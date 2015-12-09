//
//  TestTranslation.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 3/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestTranslation.h"
#include "ProcessOperatorsOMP.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorBodyVelocities.h"
#include <sstream>

void TestTranslation::_ic()
{
	Real center[3] = {.25,.25,.25};
	Real radius = .05;
	Real rhoS = 2;
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
    bool bPeriodic[3] = {false,false,false};
    
	const Real domainSize[3] = { FluidBlock::sizeX * grid->getBlocksPerDimension(0) * vInfo[0].h_gridpoint,
								 FluidBlock::sizeY * grid->getBlocksPerDimension(1) * vInfo[0].h_gridpoint,
								 FluidBlock::sizeZ * grid->getBlocksPerDimension(2) * vInfo[0].h_gridpoint};
	
	shape = new Sphere(center, radius, rhoS, (Real)2, (Real)2);
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++)
		{
			Real p[3];
			info.pos(p, ix, iy, iz);
			
			// this is for testCase==1, no effect on testCase==0
			b(ix,iy,iz).u = 1;
			b(ix,iy,iz).v = .5;
			b(ix,iy,iz).w = .25;
			
			b(ix,iy,iz).chi = shape->chi(p, info.h_gridpoint);
			
			// assume fluid with density 1
			b(ix,iy,iz).rho = shape->rho(p, info.h_gridpoint, b(ix,iy,iz).chi);
			
			// this is for testing purposes only! do it the clean way!!
			b(ix,iy,iz).p = 0;
			b(ix,iy,iz).divU = 0;
			b(ix,iy,iz).pOld = 0;
		}
	}
}

TestTranslation::TestTranslation(const int argc, const char ** argv, const int testCase, const int bpd, const double dt) : Test(argc,argv), testCase(testCase), bpd(bpd), dt(dt), nsteps(10)
{
	grid = new FluidGrid(bpd,bpd,bpd);
	
	path2file = parser("-file").asString("../data/testTranslation");
	_ic();
	
}

TestTranslation::~TestTranslation()
{
	delete grid;
}

void TestTranslation::run()
{
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	const int sizeZ = bpd * FluidBlock::sizeZ;
	
	Real u[3] = {0,0,0};
	Real lambda = 1;
	CoordinatorComputeShape coordComputeShape(shape, grid);
	CoordinatorBodyVelocities coordBodyVelocities(&u[0], &u[1], &u[2], &lambda, shape, grid);
	
	for (int step=0; step<nsteps; step++)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		if (testCase==0)
		{
			u[0] = 1;
			u[1] = .5;
			u[2] = .25;
			const Real mass = 1;
			const Real dthetadt[3] = { 0,0,0 };
			const Real J[6] = { 1,1,1,0,0,0 };
			shape->updatePosition(u, dthetadt, J, mass, dt);
		}
		else if (testCase==1)
		{
			coordBodyVelocities(dt);
		}
		else
			abort();
		
		coordComputeShape(dt);
		
		// dump
		if (step%10==0)
		{
			stringstream ss;
			ss << path2file << "-" << bpd << "-" << step << ".vti";
			
			dumper.Write(*grid, ss.str());
		}
	}
}

void TestTranslation::check()
{
	// the only thing to check here is the orientation
	Real p[3];
	shape->getCenterOfMass(p);
	cout << "Translation error X: " << p[0] - (.25 + dt*nsteps*1)   << endl;
	cout << "Translation error Y: " << p[1] - (.25 + dt*nsteps*.5)  << endl;
	cout << "Translation error Z: " << p[2] - (.25 + dt*nsteps*.25) << endl;
}