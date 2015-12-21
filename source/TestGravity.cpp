//
//  TestGravity.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 2/3/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestGravity.h"
#include "ProcessOperatorsOMP.h"
#include "CoordinatorGravity.h"
#include <sstream>
#include <cmath>

void TestGravity::_ic()
{
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
					double p[3];
					info.pos(p, ix, iy, iz);
					
					b(ix, iy, iz).rho = 0;
					b(ix, iy, iz).u   = 0;
					b(ix, iy, iz).v   = p[0];
					b(ix, iy, iz).p   = 0;
					b(ix, iy, iz).chi = 0;
					b(ix, iy, iz).tmpU = 0;
					b(ix, iy, iz).tmpV = 0;
				}
	}
	
	
#ifdef _USE_HDF_
	CoordinatorVorticity<Lab> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}

TestGravity::TestGravity(const int argc, const char ** argv, const int bpd, const double dt) : Test(argc, argv, bpd), time(0), gravity{0,-9.81,0}, dt(dt)
{
	path2file = parser("-file").asString("../data/testGravity");
	_ic();
}

TestGravity::~TestGravity()
{
}

void TestGravity::run()
{
	time = 0;
	
	const double dt = 1e-4;
	
	const int nsteps = 1000;
	
	CoordinatorGravity coordGravity(gravity, grid);
	
	for(int step=0; step<nsteps; ++step)
	{
		coordGravity(dt);
		
		time += dt;
	}
	
	
#ifdef _USE_HDF_
	CoordinatorVorticity<Lab> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << "-Final_" << bpd;
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 1, ss.str());
#endif
}

void TestGravity::check()
{
	cout << "\tErrors (uLinf, vLinf, uL1, vL1, uL2, vL2):\t";
	double uLinf = 0.;
	double vLinf = 0.;
	double vL1 = 0.;
	double uL1 = 0.;
	double uL2 = 0.;
	double vL2 = 0.;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2) reduction(max:vLinf) reduction(+:vL1) reduction(+:vL2)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy, iz);
					
					double uError = b(ix, iy, iz).u - gravity[0]*time;
					double vError = b(ix, iy, iz).v - gravity[1]*time - p[0];
					
					uLinf = max(uLinf,abs(uError));
					uL1 += abs(uError);
					uL2 += uError*uError;
					
					vLinf = max(vLinf,abs(vError));
					vL1 += abs(vError);
					vL2 += vError*vError;
				}
	}
	
	uL1 *= dh*dh*dh;
	vL1 *= dh*dh*dh;
	uL2 = sqrt(uL2*dh*dh*dh);
	vL2 = sqrt(vL2*dh*dh*dh);
	const int size = bpd * FluidBlock::sizeX;
	cout << "\t" << uLinf << "\t" << vLinf << "\t" << uL1 << "\t" << vL1 << "\t" << uL2 << "\t" << vL2 << endl;
	myfile << size << " " << uLinf << " " << vLinf << " " << uL1 << " " << vL1 << " " << uL2 << " " << vL2 << endl;
}
