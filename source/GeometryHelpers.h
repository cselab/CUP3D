//
//  GeometryHelpers.h
//  CubismUP_3D
//
//	Helper structures for WavefrontReader
//
//  Created by Christian Conti on 07/04/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_GeometryHelpers_h
#define CubismUP_3D_GeometryHelpers_h

#pragma once

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include "common.h"

using namespace std;

namespace Geometry
{
	struct Grid
	{
		vector<float> data;
		int sizeX, sizeY, sizeZ;
		int maxpoints;
		double maxsize;
		
		bool bAllocated;
		
		Grid() : bAllocated(false) {}
		
		~Grid() {}
		
		inline Grid & operator=(const Grid & g)
		{
			this->sizeX = g.sizeX;
			this->sizeY = g.sizeY;
			this->sizeZ = g.sizeZ;
			
			this->data = g.data;
			
			this->maxpoints = g.maxpoints;
			this->maxsize = g.maxsize;
			this->bAllocated = g.bAllocated;
			
			return *this;
		}
		
		inline float& operator()(int ix, int iy, int iz)
		{
			if (!bAllocated)
			{
				cout << "Geometry grid not allocated!\n";
				abort();
			}
			
			if (ix<sizeX && iy<sizeY && iz<sizeZ)
				return data[ix + iy*sizeX + iz*sizeY*sizeX];
		}
		
		inline float read(int ix, int iy, int iz, int def=1.)
		{
			if (!bAllocated)
			{
				cout << "Geometry grid not allocated!\n";
				abort();
			}
			
			if (ix>=0 && ix<sizeX &&
				iy>=0 && iy<sizeY &&
				iz>=0 && iz<sizeZ)
				return data[ix + iy*sizeX + iz*sizeY*sizeX];
			else
				return def;
		}
		
		void reset(int sx, int sy, int sz, int value=10000.)
		{
			sizeX = sx;
			sizeY = sy;
			sizeZ = sz;
			
			data.clear();
			data.resize(sizeX*sizeY*sizeZ, value);
			
			bAllocated = true;
		}
	};
	
	struct Point
	{
		double x,y,z;
		
		Point() : x(0.), y(0.), z(0.) {}
		
		Point(double px, double py, double pz)
		{
			x = px;
			y = py;
			z = pz;
		}
	};
	
	struct Quaternion
	{
		double w,x,y,z;
		
		Quaternion()
		{
			w = 1.;
			x = 0.;
			y = 0.;
			z = 0.;
		}
		
		Quaternion(double w, double x, double y, double z) : w(w), x(x), y(y), z(z) {}
		
		const double magnitude() const
		{
			return sqrt(x*x + y*y + z*z + w*w);
		}
		
		Quaternion& operator=(const Quaternion& rhs)
		{
			this->w = rhs.w;
			this->x = rhs.x;
			this->y = rhs.y;
			this->z = rhs.z;
			return *this;
		}
		
		const Quaternion operator*(const Quaternion& b)
		{
			return Quaternion(this->w * b.w - this->x * b.x - this->y * b.y - this->z * b.z,
							  this->w * b.x + this->x * b.w + this->y * b.z - this->z * b.y,
							  this->w * b.y - this->x * b.z + this->y * b.w + this->z * b.x,
							  this->w * b.z + this->x * b.y - this->y * b.x + this->z * b.w);
		}
	};
	
	struct Triangle
	{
		int a,b,c;
	};
	
	struct Properties
	{
		// quantities defined in the global system of reference
		Point com; // center of mass
		Point centroid; // centroid
		Point ut; // translational velocity
		Point dthetadt; // angular velocity (in the paper)
		double mass;
		double density;
		Point minb, maxb; // bounding box
		Quaternion q; //quaternion representing 3D rotation: qw + i qx + j qy + k qz = [qw,qx,qy,qz]
		double J[6]; // 00,11,22,01,02,12
		double rotation[3][3];
		
		Properties()
		{
			mass = 0.;
			/*
			Quaternion q1(cos(55./180.*M_PI), 0, 0, sin(55./180.*M_PI));
			Quaternion q2(cos(-12./180.*M_PI), 0,sin(-12./180.*M_PI),0);
			
			q = q1*q2;
			
			double qm = q.magnitude();
			q.w /= qm;
			q.x /= qm;
			q.y /= qm;
			q.z /= qm;
			*/
			rotation[0][0] = 1-2*(q.y*q.y + q.z*q.z);
			rotation[0][1] =   2*(q.x*q.y - q.w*q.z);
			rotation[0][2] =   2*(q.x*q.z + q.w*q.y);
			rotation[1][0] =   2*(q.x*q.y + q.w*q.z);
			rotation[1][1] = 1-2*(q.x*q.x + q.z*q.z);
			rotation[1][2] =   2*(q.y*q.z - q.w*q.x);
			rotation[2][0] =   2*(q.x*q.z - q.w*q.y);
			rotation[2][1] =   2*(q.y*q.z + q.w*q.x);
			rotation[2][2] = 1-2*(q.x*q.x + q.y*q.y);
		}
		
		// dt is the b factor of LSRK3
		void update(double dt)
		{
			// movements in global coordinates
			//cout << "Current velocity is " << ut.x << " " << ut.y << " " << ut.z << endl;
			
			// x_cm += u_T * dt
			com.x += ut.x * dt;
			com.y += ut.y * dt;
			com.z += ut.z * dt;
			
			centroid.x += ut.x * dt;
			centroid.y += ut.y * dt;
			centroid.z += ut.z * dt;
			
			minb.x += ut.x * dt;
			minb.y += ut.y * dt;
			minb.z += ut.z * dt;
			maxb.x += ut.x * dt;
			maxb.y += ut.y * dt;
			maxb.z += ut.z * dt;
			
			//cout << "Mass is " << mass << endl;
			
			// normality preserving advection (Simulation of colliding constrained rigid bodies - Kleppmann 2007 Cambridge University, p51)
			// move the correct distance on the quaternion unit ball surface, end up with normalized quaternion
			const Quaternion dq((Real)0.5*dt*(-dthetadt.x * q.x - dthetadt.y * q.y - dthetadt.z * q.z),
								(Real)0.5*dt*( dthetadt.x * q.w + dthetadt.y * q.z - dthetadt.z * q.y),
								(Real)0.5*dt*(-dthetadt.x * q.z + dthetadt.y * q.w + dthetadt.z * q.x),
								(Real)0.5*dt*( dthetadt.x * q.y - dthetadt.y * q.x + dthetadt.z * q.w));
			
			const Real dq_length = dq.magnitude();
			
			if(dq_length > numeric_limits<Real>::epsilon())
			{
				const Real tanfac = tan(dq_length)/dq_length;
				const Quaternion num(q.w + tanfac*dq.w,
									 q.x + tanfac*dq.x,
									 q.y + tanfac*dq.y,
									 q.z + tanfac*dq.z);
				
				const Real invDenum = 1./num.magnitude();
				
				q.w = num.w * invDenum;
				q.x = num.x * invDenum;
				q.y = num.y * invDenum;
				q.z = num.z * invDenum;
			}
			
			// check that the quaternion is a unit quaternion
			const Real d = q.magnitude();
			if (abs(d-1.0) > 5*numeric_limits<Real>::epsilon())
			{
				cout << "quaternion is not normalized\n";
				abort();
			}
			
			rotation[0][0] = 1-2*(q.y*q.y + q.z*q.z);
			rotation[0][1] =   2*(q.x*q.y - q.w*q.z);
			rotation[0][2] =   2*(q.x*q.z + q.w*q.y);
			rotation[1][0] =   2*(q.x*q.y + q.w*q.z);
			rotation[1][1] = 1-2*(q.x*q.x + q.z*q.z);
			rotation[1][2] =   2*(q.y*q.z - q.w*q.x);
			rotation[2][0] =   2*(q.x*q.z - q.w*q.y);
			rotation[2][1] =   2*(q.y*q.z + q.w*q.x);
			rotation[2][2] = 1-2*(q.x*q.x + q.y*q.y);
		}
	};
}

#endif
