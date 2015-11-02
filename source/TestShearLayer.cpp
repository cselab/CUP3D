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

void TestShearLayer::_getRefs(const int ix, const int iy, const int ratio, Real &u, Real &v)
{
	const int ixRef = ix * ratio;
	const int iyRef = iy * ratio;
	
	const int ilx = ixRef % FluidBlock::sizeX;
	const int ily = iyRef % FluidBlock::sizeY;
	const int bx = (ixRef-ilx) / FluidBlock::sizeX;
	const int by = (iyRef-ily) / FluidBlock::sizeY;
	
	FluidBlock& b = (*gridRef)(bx,by);
	u = b(ilx,ily).u;
	v = b(ilx,ily).v;
}

void TestShearLayer::_ic()
{
	// setup initial conditions
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
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
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
#ifndef _VARDENSITY_
				b(ix,iy).rho = 1;
#else
				const double rhoS = 7.2314;
				const double rescale = (rhoS-1)*.5;
				b(ix,iy).rho = p[1]<.5 ? 1 + rescale * (1 + tanh(r * (p[1]-.25))) : 1 + rescale * (1 + tanh(r * (.75 - p[1])));
#endif
				
				b(ix,iy).u = p[1]<.5 ? tanh(r * (p[1]-.25)) : tanh(r * (.75 - p[1]));
				b(ix,iy).v = delta * sin(2*M_PI*(p[0]+.25));
				b(ix,iy).p = 0;
				
				b(ix,iy).chi = 0;
				b(ix,iy).divU = 0;
				b(ix,iy).pOld = 0;
				
				b(ix,iy).tmpU = 0;
				b(ix,iy).tmpV = 0;
				b(ix,iy).tmp  = 0;
			}
	}
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	dumper.Write(*grid, ss.str());
	
#ifndef _THIN_
	// setup initial conditions for Ref
	vector<BlockInfo> vInfoRef = gridRef->getBlocksInfo();
	const int NRef = vInfoRef.size();
	
#pragma omp parallel for schedule(static)
	for(int i=0; i<NRef; i++)
	{
		BlockInfo infoRef = vInfoRef[i];
		FluidBlock& b = *(FluidBlock*)infoRef.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				infoRef.pos(p, ix, iy);
				
				
#ifndef _VARDENSITY_
				b(ix,iy).rho = 1;
#else
				const double rescale = (rhoS-1)*.5;
				b(ix,iy).rho = p[1]<.5 ? 1 + rescale * (1 + tanh(r * (p[1]-.25))) : 1 + rescale * (1 + tanh(r * (.75 - p[1])));
#endif
				
				b(ix,iy).u = p[1]<.5 ? tanh(r * (p[1]-.25)) : tanh(r * (.75 - p[1]));
				b(ix,iy).v = delta * sin(2*M_PI*(p[0]+.25));
				b(ix,iy).p = 0;
				
				b(ix,iy).chi = 0;
				b(ix,iy).divU = 0;
				b(ix,iy).pOld = 0;
				
				b(ix,iy).tmpU = 0;
				b(ix,iy).tmpV = 0;
				b(ix,iy).tmp  = 0;
			}
	}
	
	stringstream ssRef;
	ssRef << path2file << "-ICRef.vti";
	dumper.Write(*gridRef, ssRef.str());
#endif
}

TestShearLayer::TestShearLayer(const int argc, const char ** argv, const int bpd) : Test(argc, argv), dtCFL(0), dtFourier(0), time(0), bpd(bpd), bpdRef(32), step(0), rhoS(7.2314),
#ifndef _THIN_
nu(0.002),
#else
nu(0.0001),
#endif
endTime(1.5)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank!=0)
		omp_set_num_threads(1);
#endif // _MULTIGRID_
	
	// output settings
	path2file = parser("-file").asString("../data/testShearLayer");
	
	grid = new FluidGrid(bpd,bpd,1);
	gridRef = new FluidGrid(bpdRef,bpdRef,1);
	
	_ic();
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
#ifndef _VARDENSITY_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,1));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,rhoS));
#endif
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
#ifndef _VARDENSITY_
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid)); // need to also test with Hypre!
#else
	Real g[2] = {0,0};
	bool bSplit = false;
	const double minRho = min(1.,rhoS);
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, g, &step, bSplit, grid, rank, nprocs));
#endif
}

TestShearLayer::~TestShearLayer()
{
	delete grid;
	delete gridRef;
	
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
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = CFL*vInfo[0].h_gridpoint/abs(maxU);
			assert(!std::isnan(maxU));
			dt = min(dtCFL,dtFourier);
#ifdef _PARTICLES_
			maxA = findMaxAOMP<Lab>(vInfo,*grid);
			dtLCFL = maxA==0 ? 1e5 : LCFL/abs(maxA);
			dt = min(dt,dtLCFL);
#endif
			dt = min(dt,nextDumpTime-time);
			dt = min(dt,endTime-time);
			
			//*
			if (dt==dtFourier) cout << "Diffusion limited\n";
			else if (dt==dtCFL) cout << "Advection CFL limited\n";
			else if (dt==dtLCFL) cout << "Advection LCFL limited\n";
			else if (dt==nextDumpTime-time) cout << "dump limited\n";
			else if (dt==endTime-time) cout << "endTime limited\n";
			//*/
		}
#ifdef _MULTIGRID_
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif // _MULTIGRID_
		
		
		if (dt!=0)
		{
			for (int c=0; c<pipeline.size(); c++)
			{
#ifdef _MULTIGRID_
				MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
				if (rank == 0 || pipeline[c]->getName()=="Pressure")
					(*pipeline[c])(dt);
			}
			
			time += dt;
			step++;
			//cout << time << endl;
		}
		
		if (rank==0 && abs(time-nextDumpTime) < 10*std::numeric_limits<Real>::epsilon())
		{
			stringstream ss;
			ss << path2file << bpd << "-" << step <<  ".vti";
			dumper.Write(*grid, ss.str());
			/*
			 const int size = bpd * FluidBlock::sizeX;
			 Layer vorticity(size,size,1);
			 processOMP<Lab, OperatorDivergenceLayer>(vorticity,vInfo,*grid);
			 stringstream sVort;
			 sVort << path2file << "Vorticity-" << bpd << "-" << step << ".vti";
			 dumpLayer2VTK(step,sVort.str(),vorticity,1);
			 */
			nextDumpTime += dumpTime;
		}
		
		// check nondimensional time
		if (rank==0 && abs(time-endTime) < 10*std::numeric_limits<Real>::epsilon())
		{
			stringstream ss;
			ss << path2file << bpd << "-Final.vti";
			
			dumper.Write(*grid, ss.str());
		}
	}
	
	if (bpd==32 && rank==0)
	{
		stringstream serializedGrid;
		serializedGrid << "SerializedGrid.grid";
		DumpZBin<FluidGrid, StreamerSerialization>(*grid, serializedGrid.str(), path2file);
	}
}

void TestShearLayer::check()
{
	if (rank==0)
	{
#ifndef _THIN_
		stringstream serializedGrid;
		serializedGrid << "SerializedGrid.grid";
		ReadZBin<FluidGrid, StreamerSerialization>(*gridRef, serializedGrid.str(), path2file);
		
		cout << "\tErrors (Linf_u, L1_u, L2_u, Linf_v, L1_v, L2_v):\t";
		double Linf_u = 0.;
		double L1_u = 0.;
		double L2_u = 0.;
		double Linf_v = 0.;
		double L1_v = 0.;
		double L2_v = 0.;
		
		
		
		const int size = bpd * FluidBlock::sizeX;
		const int resRatio = bpdRef/bpd;
		
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
#pragma omp parallel for reduction(max:Linf_u) reduction(+:L1_u) reduction(+:L2_u) reduction(max:Linf_v) reduction(+:L1_v) reduction(+:L2_v)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					double p[3];
					info.pos(p, ix, iy);
					
					Real refU, refV;
					_getRefs(ix+info.index[0]*FluidBlock::sizeX,iy+info.index[1]*FluidBlock::sizeY,resRatio,refU,refV);
					
					double errorU = b(ix, iy).u - refU;
					double errorV = b(ix, iy).v - refV;
					
					Linf_u = max(Linf_u,abs(errorU));
					L1_u += abs(errorU);
					L2_u += errorU*errorU;
					
					Linf_v = max(Linf_v,abs(errorV));
					L1_v += abs(errorV);
					L2_v += errorV*errorV;
				}
		}
		
		L2_u = sqrt(L2_u)/(double)size;
		L2_v = sqrt(L2_v)/(double)size;
		L1_u /= (double)size*size;
		L1_v /= (double)size*size;
		
		stringstream ss;
		ss << path2file << "_diagnostics.dat";
		ofstream myfile(ss.str(), fstream::app);
		cout << "\t" << Linf_u << "\t" << L1_u << "\t" << L2_u << "\t" << Linf_v << "\t" << L1_v << "\t" << L2_v << endl;
		myfile << size << " " << Linf_u << " " << L1_u << " " << L2_u << " " << Linf_v << " " << L1_v << " " << L2_v << endl;
#endif
	}
}