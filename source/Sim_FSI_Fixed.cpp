//
//  Sim_FSI_Fixed.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Fixed.h"

#include "ProcessOperatorsOMP.h"
#include "OperatorDivergence.h"
#include "OperatorGradP.h"
#include "OperatorComputeShape.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorPressure.h"


void Sim_FSI_Fixed::_diagnostics()
{
	double drag = 0;
	double volume = 0;
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for schedule(static) reduction(+:drag) reduction(+:volume)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					if (b(ix,iy).chi>0)
					{
						drag += b(ix,iy,iz).u * b(ix,iy,iz).chi;
						volume += b(ix,iy,iz).chi;
					}
				}
	}
	
	drag *= dh*dh*dh*lambda;
	volume *= dh*dh*dh;
	
	const double cD = 2.*drag/(uinf*uinf*shape->getCharLength());
	
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << drag << " " << lambda << endl;
	myfile << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << dtCFL << " " << dtFourier << " " << cD << " " << lambda << endl;
}

void Sim_FSI_Fixed::_ic()
{
	CoordinatorIC coordIC(shape,uinf,grid);
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

double Sim_FSI_Fixed::_nonDimensionalTime()
{
	return 2*time*abs(uinf)/shape->getCharLength();
}

void Sim_FSI_Fixed::_outputSettings(ostream &outStream)
{
	outStream << "Fixed_FSI\n";
	outStream << "uinf " << uinf << endl;
	outStream << "re " << re << endl;
	
	Simulation_FSI::_outputSettings(outStream);
}

void Sim_FSI_Fixed::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Fixed_FSI")
	{
		cout << "Error in deserialization - Simulation_Fixed_FSI\n";
		abort();
	}
	
	// read data
	inStream >> variableName;
	assert(variableName=="uinf");
	inStream >> uinf;
	inStream >> variableName;
	assert(variableName=="re");
	inStream >> re;
	
	Simulation_FSI::_inputSettings(inStream);
}

Sim_FSI_Fixed::Sim_FSI_Fixed(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uinf(0), re(0), nu(0), dtCFL(0), dtFourier(0)
{
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a fixed body\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_FSI_Fixed::~Sim_FSI_Fixed()
{
}

void Sim_FSI_Fixed::init()
{
	Simulation_FSI::init();
	
	if (!bRestart)
	{
		// simulation settings
		uinf = parser("-uinf").asDouble(0.1);
		re = parser("-Re").asDouble(100);
		
		nu = shape->getCharLength()*uinf/re;
		
		_ic();
	}
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,1));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid));
	pipeline.push_back(new CoordinatorPenalizationFixed(shape, &lambda, grid));
	
	cout << "Coordinator/Operator ordering:\n";
	for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
}

void Sim_FSI_Fixed::simulate()
{
	const Real uBody[3] = {0,0,0};
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
	
	double nextDumpTime = time;
	double maxU = uinf;
	double maxA = 0;
	
	while (true)
	{
		// choose dt (CFL, Fourier)
		profiler.push_start("DT");
		maxU = findMaxUOMP(vInfo,*grid);
		dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu;
		dtCFL     = CFL*vInfo[0].h_gridpoint/abs(maxU);
		dt = min(dtCFL,dtFourier);
		if (dumpTime>0)
			dt = min(dt,nextDumpTime-_nonDimensionalTime());
		if (endTime>0)
			dt = min(dt,endTime-_nonDimensionalTime());
		if (verbose)
			cout << "dt (Fourier, CFL): " << dtFourier << " " << dtCFL << endl;
		profiler.pop_stop();
		
		if (dt!=0)
		{
#ifdef _DLM_
			lambda = dlm/dt;
#endif
			
			for (int c=0; c<pipeline.size(); c++)
			{
				profiler.push_start(pipeline[c]->getName());
				(*pipeline[c])(dt);
				profiler.pop_stop();
			}
			
			time += dt;
			step++;
		}
		
		
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
		
		
		if (step % 1000 == 0)
			profiler.printSummary();
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			if (verbose)
				cout << "Finished at time " << _nonDimensionalTime() << " (target end time " << endTime << ") in " << step << " step of " << nsteps << endl;
			
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
			
			exit(0);
		}
	}
}