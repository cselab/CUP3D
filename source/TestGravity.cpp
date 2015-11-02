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
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				b(ix, iy).rho = 0;
				b(ix, iy).u   = 0;
				b(ix, iy).v   = p[0];
				b(ix, iy).p   = 0;
				b(ix, iy).chi = 0;
				b(ix, iy).tmpU = 0;
				b(ix, iy).tmpV = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	
	dumper.Write(*grid, ss.str());
}

TestGravity::TestGravity(const int argc, const char ** argv, const int bpd, const double dt) : Test(argc, argv), time(0), gravity{0,-9.81}, bpd(bpd), dt(dt)
{
	grid = new FluidGrid(bpd,bpd,1);
	
	path2file = parser("-file").asString("../data/testGravity");
	_ic();
}

TestGravity::~TestGravity()
{
	delete grid;
}

void TestGravity::run()
{
	time = 0;
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	const double dt = 1e-4;
	
	const int nsteps = 1000;
	
	CoordinatorGravity coordGravity(gravity, grid);
	
	for(int step=0; step<nsteps; ++step)
	{
		coordGravity(dt);
		
		time += dt;
	}
	
	stringstream ss;
	ss << path2file << "-bpd" << bpd << ".vti";
	
	dumper.Write(*grid, ss.str());
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
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for reduction(max:uLinf) reduction(+:uL1) reduction(+:uL2) reduction(max:vLinf) reduction(+:vL1) reduction(+:vL2)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				double uError = b(ix, iy).u - gravity[0]*time;
				double vError = b(ix, iy).v - gravity[1]*time - p[0];
				
				uLinf = max(uLinf,abs(uError));
				uL1 += abs(uError);
				uL2 += uError*uError;
				
				vLinf = max(vLinf,abs(vError));
				vL1 += abs(vError);
				vL2 += vError*vError;
			}
	}
	
	uL1 *= dh*dh;
	vL1 *= dh*dh;
	uL2 = sqrt(uL2)*dh;
	vL2 = sqrt(vL2)*dh;
	const int size = bpd * FluidBlock::sizeX;
	cout << "\t" << uLinf << "\t" << vLinf << "\t" << uL1 << "\t" << vL1 << "\t" << uL2 << "\t" << vL2 << endl;
	myfile << size << " " << uLinf << " " << vLinf << " " << uL1 << " " << vL1 << " " << uL2 << " " << vL2 << endl;
}
