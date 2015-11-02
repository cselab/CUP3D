//
//  TestVarCoeffPoisson.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/23/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestVarCoeffPoisson.h"
#include <sstream>
#include <cmath>

void TestVarCoeffPoisson::_ic()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
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
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
				p[2] = p[2]*2.-1.;
				
				double x = p[0]*M_PI;
				double y = p[1]*M_PI;
				double z = p[2]*M_PI;
				b(ix, iy, iz).v   = 0;
				b(ix, iy, iz).chi = 0;
				
				// initial guess
				b(ix, iy, iz).tmp = 0;
				
				
				if (ic==0)
				{
					// constant coefficients coefficients - v1
					b(ix, iy, iz).rho = 1;
					b(ix, iy, iz).u   = 1./(4.*M_PI*M_PI)*cos(x); // expected solution
					b(ix, iy, iz).divU = -cos(x); // rhs
				}
				else if (ic==1)
				{
					// variable coefficients - v1
					b(ix, iy, iz).rho = (M_PI*(cos(x) + 2));
					b(ix, iy, iz).u   = sin(x); // expected solution - why the factor 1/4?
					b(ix, iy, iz).divU = -8*M_PI*sin(x)/((cos(x)+2)*(cos(x)+2)); // rhs
					
					/*
					b(ix, iy).rho = (M_PI*(cos(x) + 2));
					b(ix, iy).u   = cos(x); // not converging
					b(ix, iy).divU = -2*M_PI*(cos(x)+1)/((cos(x)+2)*(cos(x)+2)); // rhs
					 */
				}
				else
					abort();
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-ic" << ic << "-bpd" << bpd << "-IC.vti";
	
	dumper.Write(*grid, ss.str());
}

TestVarCoeffPoisson::TestVarCoeffPoisson(const int argc, const char ** argv, const int ic, const int bpd) : Test(argc, argv), ic(ic), bpd(bpd)
{
	grid = new FluidGrid(bpd,1,1);
	
	// output settings
	path2file = parser("-file").asString("../data/testVarCoeffPoisson");
	
	// setup initial condition
	_ic();
}

TestVarCoeffPoisson::~TestVarCoeffPoisson()
{
	delete grid;
}

void TestVarCoeffPoisson::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	stringstream ss;
	ss << path2file;
#ifdef _MULTIGRID_
	ss << "Multigrid";
	mg.setup(grid, false, 0, 1);
	mg();
#endif // _MULTIGRID_
	
 	ss << "-ic" << ic << "-bpd" << bpd << ".vti";
	
	dumper.Write(*grid, ss.str());
}

void TestVarCoeffPoisson::check()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	
	cout << "\tErrors (Linf, L1, L2):\t";
	double Linf = 0.;
	double L1 = 0.;
	double L2 = 0.;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for reduction(max:Linf) reduction(+:L1) reduction(+:L2)
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
				p[0] = p[0]*2.-1.;
				p[1] = p[1]*2.-1.;
				p[2] = p[2]*2.-1.;
				
				double x = p[0]*M_PI;
				double y = p[1]*M_PI;
				double z = p[2]*M_PI;
				
				double error = b(ix,iy,iz).tmp - b(ix,iy,iz).u;
				//if (error > 1e5)
				//	cout << error << "\t" << b(ix,iy).tmp << "\t" << b(ix,iy).u << endl;
				
				Linf = max(Linf,abs(error));
				L1 += abs(error);
				L2 += error*error;
			}
	}
	
	L1 *= dh*dh*dh;
	L2 = sqrt(L2)*dh*dh;
	cout << "\t" << Linf << "\t" << L1 << "\t" << L2 << endl;
	myfile << FluidBlock::sizeX*bpd << " " << Linf << " " << L1 << " " << L2 << endl;
}
