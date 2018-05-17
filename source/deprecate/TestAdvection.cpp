//
//  TestAdvection.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestAdvection.h"
#include "ProcessOperatorsOMP.h"
#include <sstream>
#include <cmath>

#include "CoordinatorAdvection.h"

void TestAdvection::_icBurger()
{
	const double dh = vInfo[0].h_gridpoint;
#pragma omp parallel for
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			double p[3];
			info.pos(p, ix, iy, iz);
			b(ix, iy, iz).u   = sin(p[0]*M_PI*2);
			b(ix, iy, iz).v   = 0;
			b(ix, iy, iz).w   = 0;
			b(ix, iy, iz).chi = 0;
		}
	}
	
#ifdef _USE_HDF_
	stringstream ss;
	ss << path2file << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}

TestAdvection::TestAdvection(const int argc, const char ** argv) : Test(argc, argv), time(0)
{
	dt = parser("-dt").asDouble();
	nsteps = parser("-nsteps").asInt();
	path2file = parser("-file").asString("../data/testAdvectionBurger");
	_icBurger();
}

TestAdvection::~TestAdvection()
{ }

void TestAdvection::run()
{
	int step = 0;
	CoordinatorAdvection<LabMPI> coordTransport(grid);
	while(step<nsteps) {
		coordTransport(dt);
		time += dt;
		step++;
	}
	
#ifdef _USE_HDF_
	CoordinatorVorticity<Lab> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << "-test" << testCase << "-bpd";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 1, ss.str());
#endif
}

void TestAdvection::check()
{
	const double center[3] = {.5,.5,.5};
	
	//cout << "\tErrors (uLinf, uL1, uL2):\t";
	double localLinf = 0.;
	double localL1 = 0.;
	double localL2 = 0.;
	double Linf = 0;
	double L1 = 0;
	double L2 = 0;
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	
	const double dh = vInfo[0].h_gridpoint;
	
	
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	const int sizeZ = bpd * FluidBlock::sizeZ;
	
	if (testCase==0)
	{
#pragma omp parallel for reduction(max:localLinf) reduction(+:localL1) reduction(+:localL2)
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
						
						double error;
						error = b(ix, iy).rho - sin((p[0]-time)*8.*M_PI);//*sin((p[1]+dt)*2.*M_PI);
						b(ix,iy).chi = error;
						
						localLinf = max(localLinf,abs(error));
						localL1 += abs(error);
						localL2 += error*error;
					}
		}
	}
	else
	{
#pragma omp parallel for reduction(max:localLinf) reduction(+:localL1) reduction(+:localL2)
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
						
						Real r = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
						
						double error = b(ix, iy, iz).rho - r;
						b(ix,iy).chi = error;
						
						localLinf = max(localLinf,abs(error));
						localL1 += abs(error);
						localL2 += error*error;
					}
		}
	}
	
	MPI::COMM_WORLD.Allreduce(&localLinf, &Linf, 1, MPI::DOUBLE, MPI::MAX);
	MPI::COMM_WORLD.Allreduce(&localL1, &L1, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&localL2, &L2, 1, MPI::DOUBLE, MPI::SUM);
	
	
#ifdef _USE_HDF_
	CoordinatorVorticity<Lab> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ssol;
	ssol << path2file << "-solution" << testCase << "-bpd" << bpd;
	cout << ssol.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 1, ssol.str());
#endif
	
	if (rank==0)
	{
		L1 *= dh*dh*dh;
		L2 = sqrt(L2*dh*dh*dh);
		cout << Linf << "\t" << L1 << "\t" << L2 << endl;
		myfile << sizeX << " " << Linf << " " << L1 << " " << L2 << endl;
	}
}
