//
//  Sim_FSI_Rotation.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 12/10/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Rotation.h"

#include "ProcessOperatorsOMP.h"
#include "OperatorDivergence.h"
#include "PoissonSolverScalarFFTW.h"
#include "OperatorGradP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorBodyVelocities.h"

void Sim_FSI_Rotation::_diagnostics()
{
	stringstream ss;
	ss << path2file << "_diagnostics.dat";
	ofstream myfile(ss.str(), fstream::app);
	if (verbose)
		cout << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << lambda << endl;
	myfile << step << " " << _nonDimensionalTime() << " " << bpdx << " " << dt << " " << lambda << endl;
}

void Sim_FSI_Rotation::_ic()
{
	cout << "Initial Conditions\n";
	CoordinatorIC coordIC(shape,0,grid);
	profiler.push_start(coordIC.getName());
	coordIC(0);
	
#ifdef _USE_HDF_
	cout << "IC Dump\n";
	CoordinatorVorticity<Lab> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
	profiler.pop_stop();
}

double Sim_FSI_Rotation::_nonDimensionalTime()
{
	return time;
}

void Sim_FSI_Rotation::_outputSettings(ostream &outStream)
{
	outStream << "Rotation_FSI\n";
	outStream << "uBody " << omegaBody[0] << endl;
	outStream << "vBody " << omegaBody[1] << endl;
	outStream << "wBody " << omegaBody[2] << endl;
	outStream << "re " << re << endl;
	outStream << "nu " << nu << endl;
	
	Simulation_FSI::_outputSettings(outStream);
}

void Sim_FSI_Rotation::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Rotation_FSI")
	{
		cout << "Error in deserialization - Simulation_Rotation_FSI\n";
		abort();
	}
	
	// read data
	inStream >> variableName;
	assert(variableName=="omegaBodyX");
	inStream >> omegaBody[0];
	inStream >> variableName;
	assert(variableName=="omegaBodyY");
	inStream >> omegaBody[1];
	inStream >> variableName;
	assert(variableName=="omegaBodyZ");
	inStream >> omegaBody[2];
	inStream >> variableName;
	assert(variableName=="re");
	inStream >> re;
	inStream >> variableName;
	assert(variableName=="nu");
	inStream >> nu;
	
	Simulation_FSI::_inputSettings(inStream);
}

Sim_FSI_Rotation::Sim_FSI_Rotation(const int argc, const char ** argv) : Simulation_FSI(argc, argv), omegaBody{0,0,0}, uBody{0,0,0}, re(0), nu(0), dtBody(0), dtCFL(0), dtFourier(0)
{	
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow with a rotating body\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_FSI_Rotation::~Sim_FSI_Rotation()
{
}

void Sim_FSI_Rotation::init()
{
	Simulation_FSI::init();
	
	if (!bRestart)
	{
		// simulation settings
		re = parser("-Re").asDouble(100);
		nu = parser("-nu").asDouble(1e-2);
		
		profiler.push_start("Geometry");
		delete shape;
		lambda = parser("-lambda").asDouble(1e5);
		dlm = parser("-dlm").asDouble(1.);
		
		double rhoS = parser("-rhoS").asDouble(1);
		const Real aspectRatio = 1;
		cout << "WARNING - Aspect ratio for correct positioning of sphere not implemented yet\n";
		Real center[3] = {parser("-xpos").asDouble(.5*aspectRatio),parser("-ypos").asDouble(.85),parser("-zpos").asDouble(.5*aspectRatio)};
		
		string shapeType = parser("-shape").asString("sphere");
		const int eps = 2;
		if (shapeType=="sphere")
		{
			Real radius = parser("-radius").asDouble(0.1);
			shape = new Sphere(center, radius, rhoS, eps, eps);
		}
		else if (shapeType=="samara")
		{
			const Real center[3] = {.5,.5,.5};
			const Real moll = 2;
			const int gridsize = 1024;
			const Real scale = .2;
			const Real tx = .4;
			const Real ty = .5;
			const Real tz = .4;
			Geometry::Quaternion q1(cos(.0625*M_PI), 0, 0, sin(.0625*M_PI));
			Geometry::Quaternion q2(cos(.0625*M_PI), sin(.0625*M_PI), 0, 0);
			Geometry::Quaternion q = q2*q1;//q1*q2;
			const Real isosurface = parser("-isosurface").asDouble(.004);
			const Real charSize = 0.06;
			
			const string filename = "/users/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
			shape = new GeometryMesh(filename, gridsize, isosurface, center, charSize, rhoS, moll, moll, scale, tx, ty, tz, q);
		}
		else if (shapeType=="triangle")
		{
			const Real center[3] = {.5,.5,.5};
			const Real moll = 2;
			const int gridsize = 512;
			const Real scale = .2;
			const Real tx = .4;
			const Real ty = .5;
			const Real tz = .4;
			Geometry::Quaternion q1(cos(.0625*M_PI), 0, 0, sin(.0625*M_PI));
			Geometry::Quaternion q2(cos(.0625*M_PI), sin(.0625*M_PI), 0, 0);
			Geometry::Quaternion q = q2*q1;//q1*q2;
			const Real isosurface = parser("-isosurface").asDouble(.004);
			const Real charSize = 0.06;
			
			const string filename = "/users/cconti/CubismUP_3D/launch/geometries/TriangleSub.obj";
			cout << "Loading geometry\n";
			shape = new GeometryMesh(filename, gridsize, isosurface, center, charSize, rhoS, moll, moll, scale, tx, ty, tz, q);
			cout << "Done\n";
		}
		else
		{
			cout << "Error - this shape (" << shapeType << ") is not currently implemented! Aborting now\n";
			abort();
		}
		profiler.pop_stop();
		
		_ic();
	}
	else
	{
		cout << "restart - Simulation_Rotation\n";
		abort();
	}
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(grid,1));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu, grid));
	pipeline.push_back(new CoordinatorPressureSimple<Lab>(grid));
	pipeline.push_back(new CoordinatorBodyVelocitiesForcedRot(&lambda, shape, grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &uBody[2], shape, &lambda, grid));
	pipeline.push_back(new CoordinatorComputeShape(shape, grid));
	
	cout << "Coordinator/Operator ordering:\n";
	for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
	
	assert(uBody[1] == 0);
}

void Sim_FSI_Rotation::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
	
	double nextDumpTime = time;
	double maxU = 0;
	double maxA = 0;
	
	while (true)
	{
		// choose dt (CFL, Fourier)
		profiler.push_start("DT");
		maxU = findMaxUOMP(vInfo,*grid);
		dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu*min(shape->getRhoS(),(Real)1);
		dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
		dtBody    = max(max(abs(uBody[0]),abs(uBody[1])),abs(uBody[2]))==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/max(max(abs(uBody[0]),abs(uBody[1])),abs(uBody[2]));
		dt = min(min(dtCFL,dtFourier),dtBody);
		if (dumpTime>0)
			dt = min(dt,nextDumpTime-_nonDimensionalTime());
		if (endTime>0)
			dt = min(dt,endTime-_nonDimensionalTime());
		if (verbose)
			cout << "dt (Fourier, CFL, body): " << dtFourier << " " << dtCFL << " " << dtBody << endl;
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
		
		
		if(step % 1000 == 0)
			profiler.printSummary();
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
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
			
			exit(0);
		}
	}
}