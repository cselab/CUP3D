//
//  TestPressure.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestPressure.h"
#include "OperatorDivergence.h"
#include "OperatorDiffusion.h"
#include "OperatorGradP.h"
#include "PoissonSolverScalarFFTW.h"
#include "PoissonSolverScalarFFTW_MPI.h"
#include "ProcessOperatorsOMP.h"
#include <sstream>
#include <cmath>
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif // _MULTIGRID_

class BS4
{
public:
	static inline Real eval(Real x)
	{
		const Real t = fabs(x);
		
		if (t>2) return 0;
		
		if (t>1) return pow(2-t,3)/6;
		
		return (1 + 3*(1-t)*(1 + (1-t)*(1 - (1-t))))/6;
	}
};

void TestPressure::_ic()
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
					p[0] = p[0]*2.-1.;
					p[1] = p[1]*2.-1.;
					p[2] = p[2]*2.-1.;
					
					if (ic==0)
					{
						double x = p[0]*M_PI;
						double y = p[1]*M_PI;
						double z = p[2]*M_PI;
						
						b(ix, iy, iz).u   = 1./(4.*M_PI*M_PI)*cos(x); // expected solution
						b(ix, iy, iz).divU = -cos(x); // rhs
					}
					else if (ic==1)
					{
						const Real IrI  = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
						const double strength = 100./(1+IrI*IrI)*BS4::eval(IrI/0.5*2.5);
						
						b(ix, iy, iz).rho = 1./(4.*M_PI*M_PI)*cos(p[0]*M_PI);
						b(ix, iy, iz).u   = -p[1]*strength;
						b(ix, iy, iz).v   =  p[0]*strength;
						b(ix, iy, iz).w   = 0;
						b(ix, iy, iz).chi = 0;
						
						b(ix, iy, iz).divU = -cos(p[0]*M_PI);
					}
					else if (ic==2)
					{
						// this is a test for:
						//	0-dirichlet on x=0
						//	0-neumann on x=1
						const int size = 1/dh;
						const int bx = info.index[0]*FluidBlock::sizeX;
						const int by = info.index[1]*FluidBlock::sizeY;
						const int bz = info.index[2]*FluidBlock::sizeZ;
						p[0] = (bx+ix+.5)/(double)size;
						p[1] = (by+iy+.5)/(double)size;
						p[2] = (bz+iz+.5)/(double)size;
						double x = 4*p[0]*M_PI;
						double y = 3*p[1]*M_PI_2;
						double z = 4*p[2]*M_PI;
						//b(ix,iy).divU = 81*M_PI_2*M_PI_2 * cos(y);
						//b(ix,iy).divU = -64*M_PI_2*M_PI_2 * cos(x);
						//b(ix,iy).divU = -9*M_PI_2*M_PI_2 * cos(y) + -64*M_PI_2*M_PI_2 * sin(x);
						b(ix,iy,iz).divU = -(64+64+9)*M_PI_2*M_PI_2 * cos(y) * sin(x) * sin(z);
						b(ix,iy,iz).rho = b(ix,iy,iz).divU;
						//b(ix,iy).u = -cos(y);
						//b(ix,iy).u = cos(x);
						//b(ix,iy).u = cos(y)+sin(x);
						b(ix,iy,iz).u = cos(y)*sin(x)*sin(z);
					}
				}
	}
	
#ifdef _USE_HDF_
	stringstream ss;
	ss << path2file << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}

TestPressure::TestPressure(const int argc, const char ** argv, const int solver, const int ic, const int bpd, const double dt) : Test(argc, argv, bpd), solver(solver), ic(ic), dt(dt)
{
	// output settings
	path2file = parser("-file").asString("../data/testPressure");
	
	// setup initial condition	
	_ic();
}

TestPressure::~TestPressure()
{
}

void TestPressure::run()
{
	const int size = bpd * FluidBlock::sizeX;
	
	if (ic==1 && rank==0)
		processOMP<LabMPI, OperatorDivergence>(dt, vInfo, *grid);
	
	if (solver==0)
	{
		if (ic!=2)
		{
			PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> pressureSolver(NTHREADS,*grid);
			pressureSolver.solve(*grid);
		}
		else
		{
			PoissonSolverScalarFFTW_MPI_DCT<FluidGridMPI, StreamerDiv> pressureSolver(NTHREADS,*grid);
			pressureSolver.solve(*grid);
		}
	}
	else if (solver==1)
	{
		cout << "spectral not supported - aborting now!\n";
		abort();
	}
#ifdef _MULTIGRID_
	else if (solver==2)
	{
		const bool bConstantCoefficients = true;
		MultigridHypre mg;
		mg.setup(grid, bConstantCoefficients, rank, nprocs);
		mg();
	}
	else if (solver==3)
	{
		cout << "Don't use it for this tests\n";
		abort();
		const bool bConstantCoefficients = false;
		MultigridHypre mg;
		mg.setup(grid, bConstantCoefficients, rank, nprocs);
		mg();
	}
#endif // _MULTIGRID_
	
	if (ic==1 && rank==0)
		processOMP<LabMPI, OperatorGradP>(dt, vInfo, *grid);
	
#ifdef _USE_HDF_
	stringstream ss;
	ss << path2file << "-Final";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 1, ss.str());
#endif
}

void TestPressure::check()
{
	const int size = bpd * FluidBlock::sizeX;
	
	//cout << "\tErrors (Linf, L1, L2):\t";
	double localLinf = 0.;
	double localL1 = 0.;
	double localL2 = 0.;
	double Linf = 0;
	double L1 = 0;
	double L2 = 0;
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for reduction(max:localLinf) reduction(+:localL1) reduction(+:localL2)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double error=0;
					error = b(ix,iy,iz).divU - b(ix,iy,iz).u;
					
					localLinf = max(localLinf,abs(error));
					localL1 += abs(error);
					localL2 += error*error;
				}
	}
	
	MPI::COMM_WORLD.Allreduce(&localLinf, &Linf, 1, MPI_DOUBLE, MPI_MAX);
	MPI::COMM_WORLD.Allreduce(&localL1, &L1, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&localL2, &L2, 1, MPI_DOUBLE, MPI_SUM);
	
	if (rank==0)
	{
		stringstream ss;
		ss << path2file << "_diagnostics.dat";
		ofstream myfile(ss.str(), fstream::app);
		L1 *= dh*dh*dh;
		L2 = sqrt(L2*dh*dh*dh);
		cout << "\t" << Linf << "\t" << L1 << "\t" << L2 << endl;
		myfile << size << " " << Linf << " " << L1 << " " << L2 << endl;
	}
}
