//
//  TestPoiseuille.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 5/13/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestPoiseuille.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPressure.h"
#include "CoordinatorPressureGradient.h"

void TestPoiseuille::_analytical(Real x, Real y, Real z, double t, Real &u, Real &v, Real &w, Real &p)
{
	u = 0;
	v = 0;
	w = 0;
	p = 0;
}

void TestPoiseuille::_ic()
{
	// setup initial conditions
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	const int N = vInfo.size();
	
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
				
				b(ix,iy,iz).rho = 1.;
				
				b(ix,iy,iz).u = .1;
				b(ix,iy,iz).v = 0;
				b(ix,iy,iz).w = 0;
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
	
	stringstream ss;
	ss << path2file << "-IC.vti";
	dumper.Write(*grid, ss.str());
}

TestPoiseuille::TestPoiseuille(const int argc, const char ** argv, const int bpd) : Test(argc, argv), dtCFL(0), dtFourier(0), time(0), bpd(bpd), nu(0.01), endTime(10), step(0), rank(0), nprocs(1)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#endif // _MULTIGRID_
	
	// output settings
	path2file = parser("-file").asString("../data/testPoiseuille");
	
	grid = new FluidGrid(bpd,bpd,bpd);
	
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
		_ic();
	
	Real gradient[3] = {1.,0.,0.};
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorPressureGradient(gradient,grid));
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,1));
#endif
#ifndef _MULTIGRID_
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid)); // need to also test with Hypre!
#else
	Real g[3] = {0,0,0};
	bool bSplit = false;
	pipeline.push_back(new CoordinatorPressure<Lab>(1, g, &step, bSplit, grid, rank, nprocs));
#endif
}

TestPoiseuille::~TestPoiseuille()
{
	delete grid;
	
	while(!pipeline.empty())
	{
		GenericCoordinator * g = pipeline.back();
		pipeline.pop_back();
		delete g;
	}
}

void TestPoiseuille::run()
{
	double maxU = 0;
	double maxA = 0;
	double dt = 0;
	const double CFL = .1;//25;//0.5;//
	const double LCFL = .1;
	
	while (true)
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
			if (endTime>0)
				dt = min(dt,endTime-time);
		}
#ifdef _MULTIGRID_
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif // _MULTIGRID_
		
		
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
		
		if (step%10==0)
		{
			stringstream sstmp;
			sstmp << path2file << bpd << "-" << step << ".vti";
			dumper.Write(*grid, sstmp.str());
		}
		
		// check nondimensional time
		if (rank==0 && abs(time-endTime) < 10*std::numeric_limits<Real>::epsilon())
		{
			stringstream ss;
			ss << path2file << bpd << "-Final.vti";
			
			dumper.Write(*grid, ss.str());
			
			return;
		}
	}
}

void TestPoiseuille::check()
{
	if (rank==0)
	{
#ifdef _MULTIGRID_
		cout << "Hypre - double";
#else
		cout << "FFTW - single";
#endif
		cout << "\tErrors (Linf_u, L1_u, L2_u, Linf_v, L1_v, L2_v, Linf_p, L1_p, L2_p):\t";
		double Linf_u = 0.;
		double L1_u = 0.;
		double L2_u = 0.;
		double Linf_v = 0.;
		double L1_v = 0.;
		double L2_v = 0.;
		double Linf_p = 0.;
		double L1_p = 0.;
		double L2_p = 0.;
		
		const int size = bpd * FluidBlock::sizeX;
		
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		
#pragma omp parallel for reduction(max:Linf_u) reduction(+:L1_u) reduction(+:L2_u) reduction(max:Linf_v) reduction(+:L1_v) reduction(+:L2_v) reduction(max:Linf_p) reduction(+:L1_p) reduction(+:L2_p)
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
					
					Real velU, velV, velW, pressure;
					_analytical(p[0], p[1], p[2], time, velU, velV, velW, pressure);
					
					double errorU = b(ix, iy, iz).u - velU;
					double errorV = b(ix, iy, iz).v - velV;
					double errorP = b(ix, iy, iz).p - pressure;
					
					Linf_u = max(Linf_u,abs(errorU));
					L1_u += abs(errorU);
					L2_u += errorU*errorU;
					
					Linf_v = max(Linf_v,abs(errorV));
					L1_v += abs(errorV);
					L2_v += errorV*errorV;
					
					Linf_p = max(Linf_p,abs(errorP));
					L1_p += abs(errorP);
					L2_p += errorP*errorP;
				}
		}
		
		const double invh3 = 1./((double)size*size*size);
		L2_u = sqrt(L2_u*invh3);
		L2_v = sqrt(L2_v*invh3);
		L2_p = sqrt(L2_p*invh3);
		L1_u *= invh3;
		L1_v *= invh3;
		L1_p *= invh3;
		
		stringstream ss;
		ss << path2file << "_diagnostics.dat";
		ofstream myfile(ss.str(), fstream::app);
		cout << "\t" << Linf_u << "\t" << L1_u << "\t" << L2_u << "\t" << Linf_v << "\t" << L1_v << "\t" << L2_v << "\t" << Linf_p << "\t" << L1_p << "\t" << L2_p << endl;
		myfile << size << "\t" << Linf_u << "\t" << L1_u << "\t" << L2_u << "\t" << Linf_v << "\t" << L1_v << "\t" << L2_v << "\t" << Linf_p << "\t" << L1_p << "\t" << L2_p << endl;
	}
}