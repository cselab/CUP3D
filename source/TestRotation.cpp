//
//  TestRotation.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 3/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestRotation.h"
#include "ProcessOperatorsOMP.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorBodyVelocities.h"
#include <sstream>
#include <cmath>

void TestRotation::_ic()
{
	const Real center[3] = {.5,.5,.5};
	const Real rhoS = 2;
	const Real moll = 2;
	const int gridsize = 128;
	const Real scale = .15;
	const Real tx = .425222;
	const Real ty = .455217;
	const Real tz = .377632;
	//Geometry::Quaternion q(1, 0, 0, 0);
	//Geometry::Quaternion q(cos(.125*M_PI), sin(.125*M_PI), 0, 0);
	Geometry::Quaternion q(cos(.125*M_PI), 0, sin(.125*M_PI), 0);
	const string filename = "/cluster/home/infk/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
	shape = new GeometryMesh(filename, gridsize, .0, center, rhoS, moll, moll, scale, tx, ty, tz, q);
	Real com[3];
	shape->getCenterOfMass(com);
	
	//const Real radius = .2;
	//const Real rhoS = 15./(8.*M_PI * radius*radius*radius*radius*radius);
	//shape = new Sphere(center, radius, rhoS, moll, moll);
	
	cout << "\tGeometricCoM:\t" << com[0] << " " << com[1] << " " << com[2] << endl;
	
	const double dh = vInfo[0].h_gridpoint;
	
	Real cx = 0;
	Real cy = 0;
	Real cz = 0;
	Real vol = 0;
	Real cxG = 0;
	Real cyG = 0;
	Real czG = 0;
	Real volG = 0;
	
#pragma omp parallel for reduction(+:cx) reduction(+:cy) reduction(+:cz) reduction(+:vol)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
					
					b(ix,iy,iz).chi = shape->chi(p, dh);
					b(ix,iy,iz).rho = shape->rho(p, dh, b(ix,iy,iz).chi);
					
					// this is for testing purposes only! do it the clean way!!
					b(ix,iy,iz).p = 0;
					b(ix,iy,iz).divU = 0;
					b(ix,iy,iz).pOld = 0;
					
					const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
					cx += p[0] * rhochi * dh*dh*dh;
					cy += p[1] * rhochi * dh*dh*dh;
					cz += p[2] * rhochi * dh*dh*dh;
					vol += rhochi * dh*dh*dh;
				}
	}
	
	MPI::COMM_WORLD.Allreduce(&cx, &cxG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&cy, &cyG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&cz, &czG, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&vol, &volG, 1, MPI::DOUBLE, MPI::SUM);
	
	cxG /= volG;
	cyG /= volG;
	czG /= volG;
	
	cout << "\tGridCM:\t" << cxG << " " << cyG << " " << czG << endl;
	com[0] = cxG;
	com[1] = cyG;
	com[2] = czG;
	shape->setCenterOfMass(com);
	shape->setCentroid(com);
	Real J0 = 0;
	Real J1 = 0;
	Real J2 = 0;
	Real J3 = 0;
	Real J4 = 0;
	Real J5 = 0;
	Real J0G = 0;
	Real J1G = 0;
	Real J2G = 0;
	Real J3G = 0;
	Real J4G = 0;
	Real J5G = 0;
	
#pragma omp parallel for reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5)
	for(int i=0; i<(int)vInfo.size(); i++)
	{
		BlockInfo info = vInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		const Real h = info.h_gridpoint;
		const Real h3 = h*h*h;
		
		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
				for(int ix=0; ix<FluidBlock::sizeX; ix++)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
					
					p[0] -= com[0];
					p[1] -= com[1];
					p[2] -= com[2];
					
					// this is for testCase==1, no effect on testCase==0
					//b(ix,iy,iz).u =  p[1] * 2 * M_PI;
					//b(ix,iy,iz).v = -p[0] * 2 * M_PI;
					//b(ix,iy,iz).w = 0;
					
					b(ix,iy,iz).u = -p[2] * 2 * M_PI;
					b(ix,iy,iz).v = 0;
					b(ix,iy,iz).w =  p[0] * 2 * M_PI;
					
					//b(ix,iy,iz).u = 0;
					//b(ix,iy,iz).v =  p[2] * 2 * M_PI;
					//b(ix,iy,iz).w = -p[1] * 2 * M_PI;
					
					const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
					J0 += rhochi * (p[1]*p[1] + p[2]*p[2]) * h3; //       y^2 + z^2
					J1 += rhochi * (p[0]*p[0] + p[2]*p[2]) * h3; // x^2 +     + z^2
					J2 += rhochi * (p[0]*p[0] + p[1]*p[1]) * h3; // x^2 + y^2
					J3 -= rhochi * p[0] * p[1] * h3; // xy
					J4 -= rhochi * p[0] * p[2] * h3; // xz
					J5 -= rhochi * p[1] * p[2] * h3; // yz
				}
	}
	
	MPI::COMM_WORLD.Allreduce(&J0, &J0G, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J1, &J1G, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J2, &J2G, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J3, &J3G, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J4, &J4G, 1, MPI::DOUBLE, MPI::SUM);
	MPI::COMM_WORLD.Allreduce(&J5, &J5G, 1, MPI::DOUBLE, MPI::SUM);
	
	const Real J[6] = { J0G, J1G, J2G, J3G, J4G, J5G };
	
	shape->setInertiaMatrix0(J);
	
#ifdef _USE_HDF_
	//CoordinatorVorticity<Lab> coordVorticity(grid);
	//coordVorticity(dt);
	stringstream ss;
	ss << path2file << bpd << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}

TestRotation::TestRotation(const int argc, const char ** argv, const int testCase, const int bpd, const double dt) : Test(argc,argv,bpd), testCase(testCase), dt(1./100.), nsteps(100)
{
	path2file = parser("-file").asString("../data/testRotation");
	_ic();
	
}

TestRotation::~TestRotation()
{
}

void TestRotation::run()
{
	const int sizeX = bpd * FluidBlock::sizeX;
	const int sizeY = bpd * FluidBlock::sizeY;
	const int sizeZ = bpd * FluidBlock::sizeZ;
	
	Real u[3] = {0,0,0};
	Real lambda = 1;
	Real maxU = 0;
	CoordinatorComputeShape coordComputeShape(shape, grid);
	CoordinatorBodyVelocities coordBodyVelocities(&u[0], &u[1], &u[2], &lambda, shape, &maxU, grid);
	
	Real rotation[3][3];
	shape->getOrientation(rotation);
	cout << "Orientation: " << setprecision(4) << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
	
	for (int step=0; step<nsteps; step++)
	{
		if (testCase==0)
		{
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			const Real mass = 1;
			const Real dthetadt[3] = { 0, -2*M_PI, 0 };
			const Real J[6] = { 1,1,1,0,0,0 };
			shape->updatePosition(u, dthetadt, J, mass, dt);
		}
		else if (testCase==1)
			coordBodyVelocities(dt);
		else
			abort();
		
		coordComputeShape(dt);
		
		if (rank==0)
		{
			stringstream ss;
			ss << path2file << "_diagnostics.dat";
			ofstream myfile(ss.str(), fstream::app);
			shape->getOrientation(rotation);
			myfile << step << " " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
		}
		
		//*
		if (step%5==0)
		{
			shape->getOrientation(rotation);
			cout << "Orientation (step " << step << "): " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
			
#ifdef _USE_HDF_
			CoordinatorVorticity<Lab> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << bpd << "-" << step;
			cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
		}
		//*/
	}
}

void TestRotation::check()
{
	// the only thing to check here is the orientation
	Real rotation[3][3];
	shape->getOrientation(rotation);
	//cout << "Orientation error: " << shape->getOrientation() - M_PI/4. << endl;
	cout << "Orientation error: " << 1-rotation[0][0] << " " << 0-rotation[0][1] << " " << 0-rotation[0][2] << " " << 0-rotation[1][0] << " " << 1-rotation[1][1] << " " << 0-rotation[1][2] << " " << 0-rotation[2][0] << " " << 0-rotation[2][1] << " " << 1-rotation[2][2] << " " << endl;
	//cout << "Orientation error: " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
}