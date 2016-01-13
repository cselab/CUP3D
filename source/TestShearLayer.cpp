//
//  TestShearLayer.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 4/28/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestShearLayer.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPressure.h"

#define _THIN_
#define _VARDENSITY_

void TestShearLayer::_ic()
{
	// setup initial conditions
	const int N = vInfo.size();
	
#ifndef _THIN_
	const double r = 30;
#else
	const double r = 80;
#endif
	const double delta = 0.05;
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
					
#ifndef _VARDENSITY_
					b(ix,iy,iz).rho = 1;
#else
					const double rhoS = 7.2314;
					const double rescale = (rhoS-1)*.5;
					b(ix,iy,iz).rho = p[1]<.5 ? 1 + rescale * (1 + tanh(r * (p[1]-.25))) : 1 + rescale * (1 + tanh(r * (.75 - p[1])));
#endif
					
					b(ix,iy,iz).u = p[1]<.5 ? tanh(r * (p[1]-.25)) : tanh(r * (.75 - p[1]));
					b(ix,iy,iz).v = delta * sin(2*M_PI*(p[0]+.25)) * sin(4*M_PI*(p[2]+.25));
					b(ix,iy,iz).w = p[1]<.5 ? tanh(r * (p[1]-.25)) : tanh(r * (.75 - p[1]));
					b(ix,iy,iz).p = 0;
					
					b(ix,iy,iz).chi = 0;
					b(ix,iy,iz).divU = 0;
					b(ix,iy,iz).pOld = 0;
					
					b(ix,iy,iz).tmpU = 0;
					b(ix,iy,iz).tmpV = 0;
					b(ix,iy,iz).tmpW = 0;
					b(ix,iy,iz).tmp  = 0;
				}
	}
	
#ifdef _USE_HDF_
	CoordinatorVorticity<LabMPI> coordVorticity(grid);
	coordVorticity(0);
	stringstream ss;
	ss << path2file << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
}

TestShearLayer::TestShearLayer(const int argc, const char ** argv, const int bpd) : Test(argc, argv, bpd), dtCFL(0), dtFourier(0), time(0), step(0), rhoS(7.2314),
#ifndef _THIN_
nu(0.002),
#else
nu(0.0001),
#endif
endTime(3)
{
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank!=0)
		omp_set_num_threads(1);
	
	// output settings
	path2file = parser("-file").asString("../data/testShearLayer");
	
	_ic();
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<LabMPI>(grid));
#else
#ifndef _VARDENSITY_
	pipeline.push_back(new CoordinatorAdvection<LabMPI>(grid,1));
#else
	pipeline.push_back(new CoordinatorAdvection<LabMPI>(grid,rhoS));
#endif
#endif
	pipeline.push_back(new CoordinatorDiffusion<LabMPI>(nu, grid));
#ifndef _VARDENSITY_
	pipeline.push_back(new CoordinatorPressureSimple<LabMPI>(grid)); // need to also test with Hypre!
#else
	Real g[3] = {0,0,0};
	bool bSplit = false;
	const double minRho = min(1.,rhoS);
	pipeline.push_back(new CoordinatorPressure<LabMPI>(minRho, g, &step, bSplit, grid, rank, nprocs));
#endif
}

TestShearLayer::~TestShearLayer()
{
	while(!pipeline.empty())
	{
		GenericCoordinator * g = pipeline.back();
		pipeline.pop_back();
		delete g;
	}
}

void TestShearLayer::run()
{
	double maxU = 0;
	double maxA = 0;
	double dt = 0;
	const double CFL = parser("-CFL").asDouble(.01);
	const double LCFL = .1;
	const double dumpTime = endTime/150.;
	double nextDumpTime = dumpTime;
	
	while (time<endTime)
	{
		if (rank==0)
		{
			// choose dt (CFL, Fourier)
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = CFL*vInfo[0].h_gridpoint/abs(maxU);
			assert(!std::isnan(maxU));
			dt = min(dtCFL,dtFourier);
			dt = min(dt,nextDumpTime-time);
			dt = min(dt,endTime-time);
			
			/*
			if (dt==dtFourier) cout << "Diffusion limited\n";
			else if (dt==dtCFL) cout << "Advection CFL limited\n";
			else if (dt==dtLCFL) cout << "Advection LCFL limited\n";
			else if (dt==nextDumpTime-time) cout << "dump limited\n";
			else if (dt==endTime-time) cout << "endTime limited\n";
			//*/
		}
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if (dt!=0)
		{
			for (int c=0; c<pipeline.size(); c++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				(*pipeline[c])(dt);
			}
			
			time += dt;
			step++;
			//cout << time << endl;
		}
		
		if (rank==0 && abs(time-nextDumpTime) < 10*std::numeric_limits<Real>::epsilon())
		{
#ifdef _USE_HDF_
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << bpd << "-" << step;
			cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
			
			nextDumpTime += dumpTime;
		}
		
		// check nondimensional time
		if (rank==0 && abs(time-endTime) < 10*std::numeric_limits<Real>::epsilon())
		{
#ifdef _USE_HDF_
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << "-Final";
			cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
		}
	}
}

void TestShearLayer::check()
{
}