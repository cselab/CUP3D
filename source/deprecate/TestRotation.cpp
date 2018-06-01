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
	const Real moll = 2;
#if 0
	shape = new Sphere(center, radius, rhoS, moll, moll);
#else
	const int gridsize = 1536;
	/*
	const Real scale = .3;
	const Real tx = .5;
	const Real ty = .5;
	const Real tz = .5;
	 */
	
	const Real scale = .15;//.125;
	const Real charSize = .06;
	const Real tx = .430885;
	const Real ty = .466109;
	const Real tz = .380607;
	
	//Geometry::Quaternion q(1,0,0,0);
	Geometry::Quaternion q1(cos(.125*M_PI), 0, 0, sin(.125*M_PI));
	Geometry::Quaternion q2(cos(-.3*M_PI), 0, sin(-.3*M_PI), 0);
	Geometry::Quaternion q = q1*q2;
	//const string filename = "/users/cconti/CubismUP_3D/launch/geometries/Ellipsoid.obj";
	const string filename = "/users/cconti/CubismUP_3D/launch/geometries/Samara_v3.obj";
	const Real isosurface = 0;
	shape = new GeometryMesh(filename, gridsize, isosurface, center, charSize, rhoS, moll, moll, scale, tx, ty, tz, q);
#endif
	shape->getOrientation(orientationIC);
	
	double cx = 0;
	double cy = 0;
	double cz = 0;
	double vol = 0;
	double gcx = 0;
	double gcy = 0;
	double gcz = 0;
	double gvol = 0;
	
	double J0 = 0;
	double J1 = 0;
	double J2 = 0;
	double J3 = 0;
	double J4 = 0;
	double J5 = 0;
	double J0G = 0;
	double J1G = 0;
	double J2G = 0;
	double J3G = 0;
	double J4G = 0;
	double J5G = 0;
	
	const double dh = vInfo[0].h_gridpoint;
	
#pragma omp parallel for reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5) reduction(+:cx) reduction(+:cy) reduction(+:cz) reduction(+:vol)
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
					
					b(ix,iy,iz).u = -(p[2]-.5) * 2 * M_PI;
					b(ix,iy,iz).v = 0;
					b(ix,iy,iz).w =  (p[0]-.5) * 2 * M_PI;
					
					b(ix,iy,iz).chi = shape->chi(p, dh);
					b(ix,iy,iz).rho = shape->rho(p, dh, b(ix,iy,iz).chi);
					
					b(ix,iy,iz).p = 0;
					b(ix,iy,iz).divU = 0;
					b(ix,iy,iz).pOld = 0;
					
					const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
					
					cx += p[0] * rhochi;
					cy += p[1] * rhochi;
					cz += p[2] * rhochi;
					vol += rhochi;
					
					p[0] -= center[0];
					p[1] -= center[1];
					p[2] -= center[2];
					
					J0 += rhochi * (p[1]*p[1] + p[2]*p[2]); //       y^2 + z^2
					J1 += rhochi * (p[0]*p[0] + p[2]*p[2]); // x^2 +     + z^2
					J2 += rhochi * (p[0]*p[0] + p[1]*p[1]); // x^2 + y^2
					J3 -= rhochi * p[0] * p[1]; // xy
					J4 -= rhochi * p[0] * p[2]; // xz
					J5 -= rhochi * p[1] * p[2]; // yz
				}
	}
	
	MPI::COMM_WORLD.Allreduce(&cx, &gcx, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&cy, &gcy, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&cz, &gcz, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&vol, &gvol, 1, MPI_DOUBLE, MPI_SUM);
	
	gcx /= gvol;
	gcy /= gvol;
	gcz /= gvol;
	
	cout << "Center of mass (after IC) set to " << gcx << " " << gcy << " " << gcz << endl;
	
	MPI::COMM_WORLD.Allreduce(&J0, &J0G, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&J1, &J1G, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&J2, &J2G, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&J3, &J3G, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&J4, &J4G, 1, MPI_DOUBLE, MPI_SUM);
	MPI::COMM_WORLD.Allreduce(&J5, &J5G, 1, MPI_DOUBLE, MPI_SUM);
	
	const double h3 = dh*dh*dh;
	const double J[6] = { J0G*h3, J1G*h3, J2G*h3, J3G*h3, J4G*h3, J5G*h3 };
	
	shape->setInertiaMatrix0(J);
	
#ifdef _USE_HDF_
	CoordinatorVorticity<LabMPI> coordVorticity(grid);
	coordVorticity(dt);
	stringstream ss;
	ss << path2file << bpd << "-IC";
	cout << ss.str() << endl;
	DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, 0, ss.str());
#endif
}

TestRotation::TestRotation(const int argc, const char ** argv, const int testCase, const int bpd, const double dt) : Test(argc,argv,bpd), testCase(testCase), dt(1./100.), nsteps(100), radius(.1), rhoS(1), dthetadt{0,0,0}, u{0,0,0}
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
	
	Real lambda = 1;
	double maxU = 0;
	CoordinatorComputeShape coordComputeShape(shape, grid);
	CoordinatorBodyVelocities coordBodyVelocities(&u[0], &u[1], &u[2], &lambda, shape, &maxU, grid);
	
	//Real rotation[3][3];
	//shape->getOrientation(rotation);
	//cout << "Orientation: " << setprecision(4) << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
	
	for (int step=0; step<nsteps; step++)
	{
		if (testCase==0)
		{
			u[0] = 0;
			u[1] = 0;
			u[2] = 0;
			const Real mass = rhoS * 4./3.*radius*radius*radius*M_PI;
			const double J0 = .4 * mass * radius*radius;
			dthetadt[0] = 0;
			dthetadt[1] = -2*M_PI*J0;
			dthetadt[2] = 0;
			const double J[6] = { J0,J0,J0,0,0,0 };
			shape->updatePosition(u, dthetadt, J, mass, dt);
		}
		else if (testCase==1)
		{
			coordBodyVelocities(dt);
		}
		else
		{
			cout << "Test case does not exist\n";
			abort();
		}
		
		coordComputeShape(dt);
		/*
		if (rank==0)
		{
			stringstream ss;
			ss << path2file << "_diagnostics.dat";
			ofstream myfile(ss.str(), fstream::app);
			shape->getOrientation(rotation);
			myfile << step << " " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
		}
		*/
		/*
		if (step%10==0)
		{
			//shape->getOrientation(rotation);
			//cout << "Orientation (step " << step << "): " << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
			
#ifdef _USE_HDF_
			CoordinatorVorticity<LabMPI> coordVorticity(grid);
			coordVorticity(dt);
			stringstream ss;
			ss << path2file << bpd << "-" << step;
			cout << ss.str() << endl;
			DumpHDF5_MPI<FluidGridMPI, StreamerHDF5>(*grid, step, ss.str());
#endif
		}
		//*/
	}
	
	//shape->getOrientation(rotation);
	//cout << "Orientation: " << setprecision(4) << rotation[0][0] << " " << rotation[0][1] << " " << rotation[0][2] << " " << rotation[1][0] << " " << rotation[1][1] << " " << rotation[1][2] << " " << rotation[2][0] << " " << rotation[2][1] << " " << rotation[2][2] << " " << endl;
}

void TestRotation::check()
{
	// check inertia matrix
	const Real mass = rhoS * 4./3.*radius*radius*radius*M_PI;
	const double dthetadtRef[3] = { 0, -2*M_PI, 0 };
	const double J0 = .4 * mass * radius*radius;
	//const double Jref[6] = { J0,J0,J0,0,0,0 }; // this is J in the unrotated frame for a sphere only!
	//Real J[6] = {0,0,0,0,0,0};
	Real rotation[3][3];
	
	//Real rotation[3][3];
	//shape->getOrientation(rotation);
	//shape->getInertiaMatrix(J);
	shape->getAngularVelocity(dthetadt);
	shape->getOrientation(rotation);
	
	//cout << "Inertia Matrix Error:\t" << (J[0]-Jref[0]) << " " << (J[1]-Jref[1]) << " " << (J[2]-Jref[2]) << " " << (J[3]-Jref[3]) << " " << (J[4]-Jref[4]) << " " << (J[5]-Jref[5]) << endl;
	//cout << "Angular Velocity:\t" << dthetadt[0] << " " << dthetadt[1] << " " << dthetadt[2] << endl;
	cout << "Angular Velocity Error:\t" << (dthetadt[0]-dthetadtRef[0]) << " " << (dthetadt[1]-dthetadtRef[1]) << " " << (dthetadt[2]-dthetadtRef[2]) << endl;
	cout << "Orientation Error: " << setprecision(4) << rotation[0][0]-orientationIC[0][0] << " " << rotation[0][1]-orientationIC[0][1] << " " << rotation[0][2]-orientationIC[0][2] << " " << rotation[1][0]-orientationIC[1][0] << " " << rotation[1][1]-orientationIC[1][1] << " " << rotation[1][2]-orientationIC[1][2] << " " << rotation[2][0]-orientationIC[2][0] << " " << rotation[2][1]-orientationIC[2][1] << " " << rotation[2][2]-orientationIC[2][2] << " " << endl;
}