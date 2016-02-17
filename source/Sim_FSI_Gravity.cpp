//
//  Sim_FSI_Gravity.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/26/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Sim_FSI_Gravity.h"

#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorGravity.h"
#include "CoordinatorBodyVelocities.h"

void Sim_FSI_Gravity::_diagnostics()
{
	double drag = 0;
	double volS = 0;
	double volF = 0;
	double pMin = 10;
	double pMax = 0;
	double dragG = 0;
	double volSG = 0;
	double volFG = 0;
	double pMinG = 10;
	double pMaxG = 0;
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for schedule(static) reduction(+:drag) reduction(+:volS) reduction(+:volF) reduction(max:pMax) reduction (min:pMin)
	for(int i=0; i<vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			pMin = min(pMin,(double)b(ix,iy,iz).p);
			pMax = max(pMax,(double)b(ix,iy,iz).p);
			
			if (b(ix,iy,iz).chi>0)
			{
				drag += (b(ix,iy,iz).v-uBody[1]) * b(ix,iy,iz).chi; // this depends on the direction of movement - here vertical!
				volS += b(ix,iy,iz).chi;
				volF += (1-b(ix,iy,iz).chi);
			}
			
			if (std::isnan(b(ix,iy,iz).u) ||
				std::isnan(b(ix,iy,iz).v) ||
				std::isnan(b(ix,iy,iz).rho) ||
				std::isnan(b(ix,iy,iz).chi) ||
				std::isnan(b(ix,iy,iz).p))
			{
				cout << "NaN Error - Aborting now!\n";
				abort();
			}
		}
	}
	
	drag *= dh*dh*dh*lambda;
	volS *= dh*dh*dh;
	volF *= dh*dh*dh;
	
	MPI::COMM_WORLD.Allreduce(&drag, &dragG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&volS, &volSG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&volF, &volFG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&pMin, &pMinG, 1, MPI::DOUBLE, MPI::MIN);
	MPI::COMM_WORLD.Allreduce(&pMax, &pMaxG, 1, MPI::DOUBLE, MPI::MAX);
	
	double cD = 2*dragG/(uBody[1]*uBody[1]*shape->getCharLength());
	cD = abs(uBody[1])>0 ? cD : 1e10;
	const Real Re_uBody = shape->getCharLength()*sqrt(uBody[0]*uBody[0]+uBody[1]*uBody[1]+uBody[2]*uBody[2])/nu;
	Real center[3];
	shape->getCenterOfMass(center);
	Real rotation[3][3];
	shape->getOrientation(rotation);
	
	if (rank==0)
	{
		stringstream ss;
		ss << path2file << "_diagnostics.dat";
		ofstream myfile(ss.str(), fstream::app);
		if (verbose)
		cout << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << center[0] << " " << center[1] << " " << center[2] << " " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << uBody[0] << " " << uBody[1] << " " << uBody[2] << endl;
		myfile << step << " " << time << " " << dt << " " << bpdx << " " << lambda << " " << cD << " " << Re_uBody << " " << center[0] << " " << center[1] << " " << center[2] << " " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << uBody[0] << " " << uBody[1] << " " << uBody[2] << endl;
	}
}

void Sim_FSI_Gravity::_dumpSettings(ostream& outStream)
{
	if (rank==0)
	{
		outStream << "--------------------------------------------------------------------\n";
		outStream << "Physical Settings\n";
		outStream << "\tradius\t" << shape->getCharLength()*.5 << endl;
		outStream << "\tnu\t" << nu << endl;
		outStream << "\tRhoS\t" << shape->getRhoS() << endl;
		Real center[3];
		shape->getCenterOfMass(center);
		outStream << "\tyPos\t" << center[1] << endl;
		
		outStream << "\nSimulation Settings\n";
		outStream << "\tnsteps\t" << nsteps << endl;
		outStream << "\tTend\t" << endTime << endl;
		outStream << "\tlambda\t" << lambda << endl;
#ifdef _MULTIGRID_
		outStream << "\tsplit\t" << (bSplit ? "true" : "false") << endl;
#endif
		outStream << "\tpath2file\t" << path2file << endl;
		outStream << "\tsize x\t" << nprocsx << "x" << bpdx << "x" << FluidBlock::sizeX << endl;
		outStream << "\tsize y\t" << nprocsy << "x" << bpdy << "x" << FluidBlock::sizeY << endl;
		outStream << "\tsize z\t" << nprocsz << "x" << bpdz << "x" << FluidBlock::sizeZ << endl;
#ifdef _PERIODIC_
		outStream << "\tBC\t\tperiodic\n";
#else // _PERIODIC_
		outStream << "\tBC\tmixed\n";
#endif // _PERIODIC_
#ifdef _MULTIGRID_
		outStream << "\tPoisson\tMultigrid\n";
#endif // _MULTIGRID_
#ifdef _SPLIT_
		outStream << "\tPoisson\tFFTW Split\n";
#endif // _SPLIT_
		outStream << "--------------------------------------------------------------------\n";
	}
}

void Sim_FSI_Gravity::_ic()
{
	// setup initial conditions
	//profiler.push_start("IC");
	CoordinatorIC coordIC(shape,0,grid);
	profiler.push_start(coordIC.getName());
	coordIC(0);
	profiler.pop_stop();
	
	//profiler.push_start("DumpIC");
#ifdef _USE_HDF_
	CoordinatorVorticity<LabMPI> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << "-" << std::setfill('0') << std::setw(6) << step;
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
	//profiler.pop_stop();
	
	//profiler.push_start("Diagnostics");
	_diagnostics();
	//profiler.pop_stop();
}

double Sim_FSI_Gravity::_nonDimensionalTime()
{
	return time; // how to nondimensionalize here? based on Galileo number?
}

void Sim_FSI_Gravity::_outputSettings(ostream &outStream)
{
	outStream << "Gravity_FSI\n";
	outStream << "uBody " << uBody[0] << endl;
	outStream << "vBody " << uBody[1] << endl;
	outStream << "zBody " << uBody[2] << endl;
	outStream << "re " << re << endl;
	outStream << "nu " << nu << endl;
	outStream << "minRho " << minRho << endl;
	outStream << "bSplit " << bSplit << endl;
	outStream << "maxU " << maxU << endl; // this is needed as maxU is computed only at the end of an iteration
	
	Simulation_FSI::_outputSettings(outStream);
}

void Sim_FSI_Gravity::_inputSettings(istream& inStream)
{
	string variableName;
	
	inStream >> variableName;
	if (variableName != "Gravity_FSI")
	{
		cout << "Error in deserialization - Simulation_Gravity_FSI\n";
		abort();
	}
	
	// read data
	inStream >> variableName;
	assert(variableName=="uBody");
	inStream >> uBody[0];
	inStream >> variableName;
	assert(variableName=="vBody");
	inStream >> uBody[1];
	inStream >> variableName;
	assert(variableName=="wBody");
	inStream >> uBody[2];
	inStream >> variableName;
	assert(variableName=="re");
	inStream >> re;
	inStream >> variableName;
	assert(variableName=="nu");
	inStream >> nu;
	inStream >> variableName;
	assert(variableName=="minRho");
	inStream >> minRho;
	inStream >> variableName;
	assert(variableName=="bSplit");
	inStream >> bSplit;
	
	Simulation_FSI::_inputSettings(inStream);
}

Sim_FSI_Gravity::Sim_FSI_Gravity(const int argc, const char ** argv) : Simulation_FSI(argc, argv), uBody{0,0,0}, gravity{0,-9.81,0}, dtCFL(0), dtFourier(0), dtBody(0), re(0), nu(0), minRho(0), bSplit(false)
{
	if (rank==0)
	{
		cout << "====================================================================================================================\n";
		cout << "\t\t\tFlow past a falling body\n";
		cout << "====================================================================================================================\n";
	}
}

Sim_FSI_Gravity::~Sim_FSI_Gravity()
{
}

void Sim_FSI_Gravity::init()
{
	Simulation_FSI::init();
	
	if (!bRestart)
	{
		lambda = parser("-lambda").asDouble(1e5);
		dlm = parser("-dlm").asDouble(1.);
		double rhoS = parser("-rhoS").asDouble(1);
		
		profiler.push_start("IC Geometry");
		delete shape;
		const Real aspectRatio = 1;
		if (rank==0) cout << "WARNING - Aspect ratio for correct positioning of sphere not implemented yet\n";
		Real center[3] = {parser("-xpos").asDouble(.5*aspectRatio),parser("-ypos").asDouble(.85),parser("-zpos").asDouble(.5*aspectRatio)};
		
		shapeType = parser("-shape").asString("sphere");
		const int eps = 2;
		if (shapeType=="sphere")
		{
			Real radius = parser("-radius").asDouble(0.1);
			shape = new Sphere(center, radius, rhoS, eps, eps);
		}
		else if (shapeType=="samara")
		{
#ifndef _MOVING_FRAME_
			const Real center[3] = {.5,.5,.5};
			const Real moll = 2;
			const int gridsize = 1024;
			const Real scale = .04;
			const Real tx = .12;
			const Real ty = .9;
			const Real tz = .12;
#else
			/*
			 const Real center[3] = {.5,.5,.5};
			 const Real moll = 2;
			 const int gridsize = 1024;
			 const Real scale = .04;
			 const Real tx = .12;
			 const Real ty = .9;
			 const Real tz = .12;
			 /*/
			const Real center[3] = {.5,.5,.5};
			const Real moll = 2;
			const int gridsize = 1536;
			const Real scale = .15;//.125;
			const Real tx = .45;
			const Real ty = .3;
			const Real tz = .4;
			//*/
#endif
			//Geometry::Quaternion q1(cos(.5*M_PI), 0, 0, sin(.5*M_PI));
			//Geometry::Quaternion q2(cos(0*M_PI), sin(0*M_PI), 0, 0);
			//Geometry::Quaternion q = q1*q2;
			//Geometry::Quaternion q(1,0,0,0);
			//Geometry::Quaternion q(cos(-.25*M_PI), 0, 0, sin(-.25*M_PI));
			Geometry::Quaternion q(cos(.25*M_PI), 0, 0, sin(.25*M_PI));
			const Real isosurface = parser("-isosurface").asDouble(.0);
			
			const string filename = "/users/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
			shape = new GeometryMesh(filename, gridsize, isosurface, center, rhoS, moll, moll, scale, tx, ty, tz, q);
			
			Real com[3];
			shape->getCenterOfMass(com);
			if (rank==0) cout << "Center of mass set to " << com[0] << " " << com[1] << " " << com[2] << endl;
		}
		else
		{
			cout << "Error - this shape is not currently implemented! Aborting now\n";
			abort();
		}
		profiler.pop_stop();
		
		// simulation settings
		bSplit = parser("-split").asBool(false);
		nu = parser("-nu").asDouble(1e-2);
		minRho = min((Real)1.,shape->getRhoS());
		
		stringstream ss;
		ss << path2file << "_settings.dat";
		ofstream myfile(ss.str(), fstream::app);
		_dumpSettings(cout);
		_dumpSettings(myfile);
		
		if (rank==0)
		{
#ifndef _SPLIT_
		if (bSplit)
#endif
			cout << "Using split method with constant coefficients Poisson solver\n";
#ifndef _SPLIT_
		else
			cout << "Solving full variable coefficient Poisson equation for pressure\n";
#endif
		}
		
		maxU = 0;
		
		_ic();
	}
	else
	{
		cout << "restart - Simulation_Gravity\n";
		/*
		profiler.push_start("Restart Geometry");
		delete shape;
		const Real aspectRatio = 1;
		if (rank==0) cout << "WARNING - Aspect ratio for correct positioning of sphere not implemented yet\n";
		
		const int eps = 2;
		
		// shapeType is read from the restart file
		if (shapeType=="samara")
		{
#ifndef _MOVING_FRAME_
			const Real center[3] = {.5,.5,.5};
			const Real moll = 2;
			const int gridsize = 1024;
			const Real scale = .04;
			const Real tx = .12;
			const Real ty = .9;
			const Real tz = .12;
#else
			const Real center[3] = {.5,.5,.5}; // this should come from the shape
			const Real moll = 2;
			const int gridsize = 1536;
			const Real scale = .15;//.125; // this should come from the shape
			const Real tx = .45; // this should come from the shape
			const Real ty = .3; // this should come from the shape
			const Real tz = .4; // this should come from the shape
#endif
			//Geometry::Quaternion q1(cos(.5*M_PI), 0, 0, sin(.5*M_PI));
			//Geometry::Quaternion q2(cos(0*M_PI), sin(0*M_PI), 0, 0);
			//Geometry::Quaternion q = q1*q2;
			//Geometry::Quaternion q(1,0,0,0);
			//Geometry::Quaternion q(cos(-.25*M_PI), 0, 0, sin(-.25*M_PI));
			Geometry::Quaternion q(cos(.25*M_PI), 0, 0, sin(.25*M_PI));
			const Real isosurface = parser("-isosurface").asDouble(.0); // this should come from the shape
			
			const string filename = "/users/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
			shape = new GeometryMesh(filename, gridsize, isosurface, center, rhoS, moll, moll, scale, tx, ty, tz, q);
			
			Real com[3];
			shape->getCenterOfMass(com);
			cout << "Center of mass set to " << com[0] << " " << com[1] << " " << com[2] << endl;
		}
		else
		{
			cout << "Error - this shape is not currently implemented for restart! Aborting now\n";
			abort();
		}
		profiler.pop_stop();
		 */
	}
	
	pipeline.clear();
#ifndef _MULTIPHASE_
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0],&uBody[1],&uBody[2],grid));
#else
	pipeline.push_back(new CoordinatorAdvection<Lab>(&uBody[0],&uBody[1],&uBody[2],grid,1));
#endif
	pipeline.push_back(new CoordinatorDiffusion<Lab>(nu,&uBody[0],&uBody[1],&uBody[2],grid));
	pipeline.push_back(new CoordinatorGravity(gravity, grid));
	pipeline.push_back(new CoordinatorPressure<Lab>(minRho, gravity, &uBody[0], &uBody[1], &uBody[2], &step, bSplit, grid, rank, nprocs));
	pipeline.push_back(new CoordinatorBodyVelocities(&uBody[0], &uBody[1], &uBody[2], &lambda, shape, &maxU, grid));
	pipeline.push_back(new CoordinatorPenalization(&uBody[0], &uBody[1], &uBody[2], shape, &lambda, grid)); // also computes new shape
	//pipeline.push_back(new CoordinatorComputeShape(shape, grid));
	
	if (rank==0)
	{
		cout << "Coordinator/Operator ordering:\n";
		for (int c=0; c<pipeline.size(); c++)
		cout << "\t" << pipeline[c]->getName() << endl;
	}
}

void Sim_FSI_Gravity::simulate()
{
	const int sizeX = bpdx * FluidBlock::sizeX;
	const int sizeY = bpdy * FluidBlock::sizeY;
	const int sizeZ = bpdz * FluidBlock::sizeZ;
	
	double vOld = 0;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double nextDumpTime = time;
	
	while (true)
	{
		if (rank==0)
		{
			// choose dt (CFL, Fourier)
			//profiler.push_start("DT");
			//maxU = findMaxUOMP(vInfo,*grid);
			dtFourier = CFL*vInfo[0].h_gridpoint*vInfo[0].h_gridpoint/nu*min(shape->getRhoS(),(Real)1);
			dtCFL     = maxU==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/abs(maxU);
			dtBody    = max(max(abs(uBody[0]),abs(uBody[1])),abs(uBody[2]))==0 ? 1e5 : CFL*vInfo[0].h_gridpoint/max(max(abs(uBody[0]),abs(uBody[1])),abs(uBody[2]));
			assert(!std::isnan(maxU));
			assert(!std::isnan(uBody[0]));
			assert(!std::isnan(uBody[1]));
			assert(!std::isnan(uBody[2]));
			dt = min(min(dtCFL,dtFourier),dtBody);
			if (endTime>0)
				dt = min(dt,endTime-_nonDimensionalTime());
			if (verbose)
				cout << "dt (Fourier, CFL, body): " << dt << " " << dtFourier << " " << dtCFL << " " << dtBody << endl;
			//profiler.pop_stop();
		}
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		if (time>nextDumpTime)
		{
			bDump = true;
			//nextDumpTime += dumpTime; // this is already done in _dump
		}
		
		if (dt!=0)
		{
#ifdef _DLM_
			lambda = dlm/dt;
#endif
			for (int c=0; c<pipeline.size(); c++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if (verbose) cout << pipeline[c]->getName() << endl;
				profiler.push_start(pipeline[c]->getName());
#ifdef _MULTIGRID_
				if (rank == 0 || pipeline[c]->getName()=="Pressure")
#endif
				(*pipeline[c])(dt);
				profiler.pop_stop();
			}
			
			time += dt;
			step++;
		}
		
		
		if (rank==0)
		{
			//if (step<100)
			{
				// this still needs to be corrected to the frame of reference!
				double accM = (uBody[1]-vOld)/dt;
				vOld = uBody[1];
				double accT = (shape->getRhoS()-1)/(shape->getRhoS()+.5) * gravity[1];
				double accN = (shape->getRhoS()-1)/(shape->getRhoS()   ) * gravity[1];
				if (verbose) cout << "Acceleration with added mass (measured, expected, no added mass)\t" << accM << "\t" << accT << "\t" << accN << " " << uBody[1] << " " << vOld << endl;
				stringstream ss;
				ss << path2file << "_addedmass.dat";
				ofstream myfile(ss.str(), fstream::app);
				myfile << step << " " << time << " " << accM << " " << accT << " " << accN << endl;
			}
		}
		
		// compute diagnostics
		if (step % 1 == 0)
		{
			profiler.push_start("Diagnostics");
			_diagnostics();
			profiler.pop_stop();
		}
		
		//dump some time steps every now and then
		profiler.push_start("Dump");
		_dump(nextDumpTime);
		profiler.pop_stop();
		
		if (rank==0 && step % 100 == 0)
		{
			profiler.printSummary();
		}
		
		// check nondimensional time
		if ((endTime>0 && abs(_nonDimensionalTime()-endTime) < 10*std::numeric_limits<Real>::epsilon()) || (nsteps!=0 && step>=nsteps))
		{
			profiler.push_start("Dump");
#ifdef _USE_HDF_
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << "-" << std::setfill('0') << std::setw(6) << step;
			cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
			profiler.pop_stop();
			
			if (rank==0)
				profiler.printSummary();
			
			MPI_Finalize();
			exit(0);
		}
		
		if (verbose) cout << "Barrier endstep " << rank << " step " << step << endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
}