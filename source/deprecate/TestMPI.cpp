//
//  TestMPI.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 12/18/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#include "TestMPI.h"
#include "Definitions.h"
#include <sstream>
#include <cmath>

#include "ProcessOperatorsOMP.h"
#include "OperatorTestDerivativesMPI.h"

void TestMPI::_ic()
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
					p[0] = p[0]-.5;
					p[1] = p[1]-.5;
					p[2] = p[2]-.5;
					
					// this is a test for:
					//	0-dirichlet on x=0
					//	0-neumann on x=1
					b(ix,iy,iz).u = p[0];
					b(ix,iy,iz).v = p[1]*p[1];
					b(ix,iy,iz).w = p[2]*p[2]*p[2];
					b(ix,iy,iz).chi = p[0];
					b(ix,iy,iz).p   = p[1]*p[1];
					b(ix,iy,iz).tmp = p[2]*p[2]*p[2];
				}
	}
}

TestMPI::TestMPI(const int argc, const char ** argv, const int bpd) : Test(argc, argv, bpd)
{
	// output settings
	path2file = parser("-file").asString("../data/testMPI");
	
	cout << "Running MPI test with\n";
	cout << "\tProcess grid: " << nprocsx << "x" << nprocsy << "x" << nprocsz <<endl;
	cout << "\tLocal grid: " << bpd/nprocsx << "x" << bpd/nprocsy << "x" << bpd/nprocsz <<endl;
	
	// setup initial condition
	_ic();
}

TestMPI::~TestMPI()
{
}

void TestMPI::run()
{
	processOMP<LabMPI,OperatorTestDerivativesMPI>(0,vInfo,*grid);
}

void TestMPI::check()
{
	double localLinfU = 0.;
	double localLinfV = 0.;
	double localLinfW = 0.;
	double LinfU = 0;
	double LinfV = 0;
	double LinfW = 0;
	
#pragma omp parallel for reduction(max:localLinfU) reduction(max:localLinfV) reduction(max:localLinfW)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					const int gx = ix + info.index[0]*FluidBlock::sizeX;
					const int gy = iy + info.index[1]*FluidBlock::sizeY;
					const int gz = iz + info.index[2]*FluidBlock::sizeZ;
					
					if (gx>0 && gy>0 && gz>0 && gx<bpd*FluidBlock::sizeX-1 && gy<bpd*FluidBlock::sizeY-1 && gz<bpd*FluidBlock::sizeZ-1)
					{
						double p[3];
						info.pos(p, ix, iy, iz);
						p[0] = p[0]-.5;
						p[1] = p[1]-.5;
						p[2] = p[2]-.5;
						
						double errorU = b(ix, iy, iz).chi;
						double errorV = b(ix, iy, iz).p  -2;
						double errorW = b(ix, iy, iz).tmp-6*p[2];
						
						localLinfU = max(localLinfU,abs(errorU));
						localLinfV = max(localLinfV,abs(errorV));
						localLinfW = max(localLinfW,abs(errorW));
					}
					else
					{
						localLinfU = 0;
						localLinfV = 0;
						localLinfW = 0;
					}
				}
	}
	
	MPI::COMM_WORLD.Allreduce(&localLinfU, &LinfU, 1, MPI_DOUBLE, MPI_MAX);
	MPI::COMM_WORLD.Allreduce(&localLinfV, &LinfV, 1, MPI_DOUBLE, MPI_MAX);
	MPI::COMM_WORLD.Allreduce(&localLinfW, &LinfW, 1, MPI_DOUBLE, MPI_MAX);
	
	if (rank==0)
		cout << LinfU << "\t" << LinfV << "\t" << LinfW << endl;
	
#ifdef _USE_HDF_
	stringstream ss;
	ss << path2file << "-Final";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}