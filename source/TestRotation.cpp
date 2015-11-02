//
//  TestRotation.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 3/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestRotation.h"
#include "ProcessOperatorsOMP.h"
#include "OperatorVorticity.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorBodyVelocities.h"
#include <sstream>
#include <cmath>

void TestRotation::_ic()
{
	Real center[3] = {.5,.5,.5};
	Real semiAxis[3] = {.05,.4,.2};
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
    bool bPeriodic[3] = {false,false,false};
    
	const Real domainSize[3] = { FluidBlock::sizeX * grid->getBlocksPerDimension(0) * vInfo[0].h_gridpoint,
								 FluidBlock::sizeY * grid->getBlocksPerDimension(1) * vInfo[0].h_gridpoint,
								 FluidBlock::sizeZ * grid->getBlocksPerDimension(2) * vInfo[0].h_gridpoint};
	
	const Real moll = 2;
	const Real rho = 2;
	const Real angle = M_PI/4.;
	shape = new Ellipsoid(center, semiAxis, angle, rho, moll, moll, bPeriodic, domainSize);

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
			
			const Real dist[3] = {p[0]-center[0],p[1]-center[1],p[2]-center[2]};
			
			// this is for testCase==1, no effect on testCase==0
			b(ix,iy,iz).u = -dist[1] * M_PI / 50;
			b(ix,iy,iz).v =  dist[0] * M_PI / 50;
			b(ix,iy,iz).w = 0;
			
			b(ix,iy,iz).chi = shape->chi(p, info.h_gridpoint);
			
			// assume fluid with density 1
			b(ix,iy,iz).rho = shape->rho(p, info.h_gridpoint);
			
			// this is for testing purposes only! do it the clean way!!
			b(ix,iy,iz).p = 0;
			b(ix,iy,iz).divU = 0;
			b(ix,iy,iz).pOld = 0;
		}
	}
}

TestRotation::TestRotation(const int argc, const char ** argv, const int testCase, const int bpd, const double dt) : Test(argc,argv), testCase(testCase), bpd(bpd), dt(dt)
{
	grid = new FluidGrid(bpd,bpd,bpd);
	
	path2file = parser("-file").asString("../data/testRotation");
	_ic();
	
}

TestRotation::~TestRotation()
{
	delete grid;
}

void TestRotation::run()
{
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	const int sizeZ = bpd * FluidBlock::sizeZ;
	
	Real u[3] = {0,0,0};
	Real omega = 0;
	Real lambda = 1;
	CoordinatorComputeShape coordComputeShape(&u[0], &u[1], &u[2], &omega, shape, grid);
	CoordinatorBodyVelocities coordBodyVelocities(&u[0], &u[1], &u[2], &omega, &lambda, shape->getRhoS(), grid);
	
	for (int step=0; step<100; step++)
	{
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
		if (testCase==0)
			omega = M_PI/50.; // 1 complete turn
		else if (testCase==1)
			coordBodyVelocities(dt);
		else
			abort();
		
		coordComputeShape(dt);
		
		// dump
		if (step%5==0)
		{
			stringstream ss;
			ss << path2file << "-" << bpd << "-" << step << ".vti";
			
			dumper.Write(*grid, ss.str());
			
			Layer vorticity(sizeX,sizeY,sizeZ);
			processOMP<Lab, OperatorVorticity>(vorticity,vInfo,*grid);
			stringstream sVort;
			sVort << path2file << "Vorticity-" << bpd << "-" << step << ".vti";
			dumpLayer2VTK(step,sVort.str(),vorticity,1);
		}
	}
}

void TestRotation::check()
{
	// the only thing to check here is the orientation
	cout << "Orientation error: " << shape->getOrientation() - M_PI/4. << endl;
}