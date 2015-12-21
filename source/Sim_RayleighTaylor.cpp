//
//  Sim_RayleighTaylor.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_RayleighTaylor.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"

void Sim_RayleighTaylor::_diagnostics()
{
}

void Sim_RayleighTaylor::_ic()
{
	if (rank==0)
	{
		// setup initial conditions
		CoordinatorIC_RT coordIC(rhoS, grid);
		profiler.push_start(coordIC.getName());
		coordIC(0);
		
#ifdef _USE_HDF_
		CoordinatorVorticity<Lab> coordVorticity(grid);
		coordVorticity(dt);
		stringstream ss;
		ss << path2file << "-IC";
		cout << ss.str() << endl;
		DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
		profiler.pop_stop();
	}
}

double Sim_RayleighTaylor::_nonDimensionalTime()
{
	return time*2.;
}

void Sim_RayleighTaylor::_outputSettings(ostream &outStream)
{
	outStream << "RTI\n";
	
	Simulation_MP::_outputSettings(outStream);
}

void Sim_RayleighTaylor::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "RTI")
	{
		cout << "Error in deserialization - Simulation_RayleighTaylor\n";
		abort();
	}
	
	Simulation_MP::_inputSettings(inStream);
}

Sim_RayleighTaylor::Sim_RayleighTaylor(const int argc, const char ** argv) : Simulation_MP(argc, argv)
{
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	
	if (rank!=0)
		omp_set_num_threads(1);
	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tRayleigh-Taylor Instability\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_RayleighTaylor::~Sim_RayleighTaylor()
{
}

void Sim_RayleighTaylor::init()
{
	Simulation_MP::init();
	
	if (!bRestart)
		_ic();
	
	pipeline.clear();
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,rhoS));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, gravity, &step, bSplit, grid, rank, nprocs));
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
			cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_RayleighTaylor::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double nextDumpTime = time;
	double maxU = 0;
	double maxA = 0;
	
	while (true)
	{
		if (rank==0)
		{
			profiler.push_start("DT");
			/*
			// choose dt (CFL, Fourier)
			maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			assert(!std::isnan(maxU));
			dt = min(dtCFL,dtFourier);
			 */
			dt = 0.00025/(double)bpdx;
			 
			if (dumpTime>0)
				dt = min(dt,nextDumpTime-_nonDimensionalTime());
			if (endTime>0)
				dt = min(dt,endTime-_nonDimensionalTime());
			
			if (verbose)
				cout << "t, dt (Fourier, CFL, LCFL): " << _nonDimensionalTime() << "\t" << dt << " " << dtFourier << " " << dtCFL << " " << dtLCFL << endl;
			profiler.pop_stop();
		}
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if (dt!=0)
		{
			for (int c=0; c<pipeline.size(); c++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				
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
			// compute diagnostics
			if (step % 10 == 0)
			{
				profiler.push_start("Diagnostics");
				_diagnostics();
				profiler.pop_stop();
			}
			
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
#ifdef _USE_HDF_
				CoordinatorVorticity<Lab> coordVorticity(grid);
				coordVorticity(dt);
				stringstream ss;
				ss << path2file << "-Final";
				cout << ss.str() << endl;
				DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
				profiler.pop_stop();
				profiler.printSummary();
			}
			exit(0);
		}
	}
}