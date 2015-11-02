//
//  TestDiffusion.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestDiffusion.h"
#include "CoordinatorDiffusion.h"
#include <sstream>
#include <cmath>

double TestDiffusion::_analytical(double px, double py, double t)
{
    return sin(px*2.*freq*M_PI) * sin(py*2.*freq*M_PI) * exp(-4.*2*freq*freq*nu*M_PI*M_PI*t);
}

void TestDiffusion::_ic()
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
				
				b(ix, iy).rho = 1;
				b(ix, iy).u   = _analytical(p[0],p[1],0);
				b(ix, iy).v   = 0;
				b(ix, iy).chi = 0;
				
				b(ix, iy).tmpU = 0;
				b(ix, iy).tmpV = 0;
			}
	}
	
	
	stringstream ss;
	ss << path2file << "-IC.vti" ;
	//cout << ss.str() << endl;
	
	dumper.Write(*grid, ss.str());
}

TestDiffusion::TestDiffusion(const int argc, const char ** argv, const int bpd, const double dt, const int nsteps, const double freq) : Test(argc, argv), bpd(bpd), time(0), dt(dt), nsteps(nsteps), freq(freq)
{
	// test settings
	nu = parser("-nu").asDouble(1);
	
	// output settings
	path2file = parser("-file").asString("../data/testDiffusion");
	
	grid = new FluidGrid(bpd,bpd,1);
	
	// setup initial condition
	_ic();
}

TestDiffusion::~TestDiffusion()
{
	delete grid;
}

void TestDiffusion::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	//cout << "Using dt " << dt << " (Fourier time step: " << vInfo[0].h_gridpoint*vInfo[0].h_gridpoint*.25/nu << ")\n";
	int step = 0;
	
    if (nsteps==1)
    {
        CoordinatorDiffusion<Lab> coordDiffusion(nu, grid);
        coordDiffusion(dt);
        time += dt;
        step++;
    }
	else
    {
        CoordinatorDiffusionTimeTest coordDiffusion(nu, freq, grid);
        //CoordinatorDiffusion<Lab> coordDiffusion(nu, grid);
        
        while(step<nsteps)
        {
            coordDiffusion(dt);
            
            time += dt;
            step++;
        }
    }
	
	stringstream ss;
	ss << path2file << "-bpd" << bpd << ".vti";
	
	dumper.Write(*grid, ss.str());
}

void TestDiffusion::check()
{
	//cout << "\tErrors (Linf, L1, L2):\t";
	double Linf = 0.;
	double L1 = 0.;
	double L2 = 0.;
	
	const int size = bpd * FluidBlock::sizeX;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
#pragma omp parallel for reduction(max:Linf) reduction(+:L1) reduction(+:L2)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++)
			{
				double p[3];
				info.pos(p, ix, iy);
				
				double error = b(ix, iy).u - _analytical(p[0],p[1],time);
				Linf = max(Linf,abs(error));
				L1 += abs(error);
				L2 += error*error;
			}
	}
	
	L2 = sqrt(L2)/(double)size;
	L1 /= (double)size*size;
	cout << Linf << "\t" << L1 << "\t" << L2 << endl;
	myfile << size << " " << dt << " " << Linf << " " << L1 << " " << L2 << endl;
}