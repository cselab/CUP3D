//
//  Sim_Jet.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_Jet.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"

void Sim_Jet::_diagnostics()
{
}

void Sim_Jet::_ic()
{
	if (rank==0)
	{
		profiler.push_start("IC - Jet");
		
		// setup initial conditions
		vector<BlockInfo> vInfo = grid->getBlocksInfo();
		const int N = vInfo.size();
		
		const double width = 0.01;
		const double width2 = 0.1;
		const double ampl = 0.05;
		
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
					
					if (p[1]<.5-.5*width2 || p[1]>.5+.5*width2)
						b(ix,iy,iz).rho = 1.;
					else
						b(ix,iy,iz).rho = rhoS;
					
					double aux = (width2 - 2. * abs(p[1]-.5)) / (4. * width);
					b(ix,iy,iz).u = .5 * (1 + tanh(aux)) * (1 + ampl * sin(8. * M_PI * p[0]));
					b(ix,iy,iz).v = 0;
					b(ix,iy,iz).w = 0;
					b(ix,iy,iz).chi = 0;
					
					b(ix,iy,iz).p = 0;
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
		profiler.pop_stop();
	}
}

double Sim_Jet::_nonDimensionalTime()
{
	return time; // how to nondimensionalize here? based on Galileo number?
}

void Sim_Jet::_outputSettings(ostream &outStream)
{
	outStream << "jet\n";
	
	Simulation_MP::_outputSettings(outStream);
}

void Sim_Jet::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "jet")
	{
		cout << "Error in deserialization - Simulation_Jet\n";
		abort();
	}
	
	Simulation_MP::_inputSettings(inStream);
}

Sim_Jet::Sim_Jet(const int argc, const char ** argv) : Simulation_MP(argc, argv)
{
#ifdef _MULTIGRID_
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank!=0)
		omp_set_num_threads(1);
#endif // _MULTIGRID_
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tJet\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_Jet::~Sim_Jet()
{
}

void Sim_Jet::init()
{
	Simulation_MP::init();
	
	nu = 0.00025;
	
	_ic();
	
	Real g[2] = {0,0};
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,rhoS));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, g, &step, bSplit, grid, rank, nprocs));
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
			cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_Jet::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
	
	double vOld = 0;
	
#ifdef _MULTIGRID_
	MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
	double nextDumpTime = time;
	double maxU = 0;
	double maxA = 0;
	
	while (true)
	{
		if (rank==0)
		{
			vector<BlockInfo> vInfo = grid->getBlocksInfo();
			
			// choose dt (CFL, Fourier)
			profiler.push_start("DT");
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			assert(!std::isnan(maxU));
			dt = min(dtCFL,dtFourier);
#ifdef _PARTICLES_
			maxA = findMaxAOMP<Lab>(vInfo,*grid);
			dtLCFL = maxA==0 ? 1e5 : LCFL/abs(maxA);
			dt = min(dt,dtLCFL);
#endif
			if (dumpTime>0)
				dt = min(dt,nextDumpTime-_nonDimensionalTime());
			if (endTime>0)
				dt = min(dt,endTime-_nonDimensionalTime());
			if (verbose)
				cout << "dt (Fourier, CFL): " << dt << " " << dtFourier << " " << dtCFL << endl;
			profiler.pop_stop();
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
				profiler.push_start(pipeline[c]->getName());
				if (rank == 0 || pipeline[c]->getName()=="Pressure")
					(*pipeline[c])(dt);
				profiler.pop_stop();
			}
			
			time += dt;
			step++;
		}
		
		if (rank==0)
		{
			/*
			// compute diagnostics
			if (step % 10 == 0)
			{
				profiler.push_start("Diagnostics");
				_diagnostics();
				profiler.pop_stop();
			}
			 */
			
			//dump some time steps every now and then
			profiler.push_start("Dump");
			_dump(nextDumpTime);
			profiler.pop_stop();
			
			if (step % 100 == 0)
				profiler.printSummary();
		}
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			if (rank==0)
			{
				profiler.push_start("Dump");
				stringstream ss;
				ss << path2file << "-Final.vti";
				cout << ss.str() << endl;
				
				dumper.Write(*grid, ss.str());
				profiler.pop_stop();
				
				profiler.printSummary();
			}
			exit(0);
		}
	}
}