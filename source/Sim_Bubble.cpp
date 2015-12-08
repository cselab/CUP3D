//
//  Sim_Bubble.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 4/10/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_Bubble.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"

void Sim_Bubble::_diagnostics()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	double vBubble = 0;
	double volume = 0;
	
#pragma omp parallel for schedule(static) reduction(+:vBubble) reduction(+:volume)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				if (abs(b(ix,iy,iz).rho-rhoS) < 10*std::numeric_limits<Real>::epsilon())
				{
					vBubble += b(ix,iy,iz).v;
					volume += 1;
				}
				
				if (std::isnan(b(ix,iy,iz).u) ||
					std::isnan(b(ix,iy,iz).v) ||
					std::isnan(b(ix,iy,iz).w) ||
					std::isnan(b(ix,iy,iz).rho) ||
					std::isnan(b(ix,iy,iz).chi) ||
					std::isnan(b(ix,iy,iz).p))
				{
					cout << "NaN Error - Aborting now!\n";
					abort();
				}
			}
	}
	
	vBubble /= volume;
	volume *= vInfo[0].h_gridpoint*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint;
	double vPredicted = time*gravity[1]*(rhoS-1)/(rhoS+.5);
	double vPredictedNoAM = time*gravity[1]*(rhoS-1)/(rhoS);
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << time << " " << vBubble << " " << vPredicted << " " << vPredictedNoAM << " " << volume << endl;
	myfile << step << " " << time << " " << vBubble << " " << vPredicted << " " << vPredictedNoAM << " " << volume << endl;
}

void Sim_Bubble::_ic()
{
#ifdef _MULTIGRID_
	if (rank==0)
#endif // _MULTIGRID_
	{
		// setup initial conditions
		Shape * shape;
		Real radius = parser("-radius").asDouble(0.1);
		Real centerOfMass[3];
		centerOfMass[0] = .5;
		centerOfMass[2] = .5;
		if (rhoS >= 1)
			centerOfMass[1] = .85;
		else
			centerOfMass[1] = .15;
		
		shape = new Sphere(centerOfMass, radius, rhoS, 2, 2);
		CoordinatorIC coordIC(shape,0,grid);
		profiler.push_start(coordIC.getName());
		coordIC(0);
		
		stringstream ss;
		ss << path2file << "-IC.vti";
		dumper.Write(*grid, ss.str());
		profiler.pop_stop();
	}
}

double Sim_Bubble::_nonDimensionalTime()
{
	return time; // how to nondimensionalize here? based on Galileo number?
}

void Sim_Bubble::_outputSettings(ostream &outStream)
{
	outStream << "Bubble\n";
	
	Simulation_MP::_outputSettings(outStream);
}

void Sim_Bubble::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Bubble")
	{
		cout << "Error in deserialization - Simulation_Bubble\n";
		abort();
	}
	
	Simulation_MP::_inputSettings(inStream);
}

Sim_Bubble::Sim_Bubble(const int argc, const char ** argv) : Simulation_MP(argc, argv)
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
		cout << "\t\t\tBubble flow\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_Bubble::~Sim_Bubble()
{
}

void Sim_Bubble::init()
{
	Simulation_MP::init();
	
	_ic();
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,rhoS));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, gravity, &step, bSplit, grid, rank, nprocs)); // this should only compute the dynamic pressure
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
			cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_Bubble::simulate()
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
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu*min((Real)rhoS,(Real)1);
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