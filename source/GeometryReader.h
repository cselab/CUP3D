//
//  GeometryReader.h
//  CubismUP_3D
//
//	Abstract class for geometry reader
//
//  Created by Christian Conti on 02/09/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__GeometryReader__
#define __CubismUP_3D__GeometryReader__

#pragma once

#include <fstream>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <cassert>

#include "GeometryHelpers.h"

class GeometryReader
{
protected:
	// these can be put into an enum type to simplify
	bool bGeometryLoaded;
	bool bSDFComputed;
	
	string filename;
	
	const Real baseScale;
	const Geometry::Point baseTranslation;
	
	bool bVerbose;
	
	Geometry::Properties & properties;
	int gridsize;
	
	virtual void parse(string filename) = 0;
	
public:
	// this is the grid on which to interpolate the mesh,
	// it could be eventually placed outside of this class and
	// be a cubism grid
	Geometry::Grid grid;
	
	// constructor loading OBJ file
	GeometryReader(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor);
	
	// clean up
	~GeometryReader();
	
	virtual void load(const string filename) = 0;
	
	virtual void sdf(double comx, double comy, double comz) = 0;
	
	virtual void serialize(string path) const;
	virtual void deserialize(string path);
	
	// check if a point lies inside the body
	virtual int isPointInside(int ix, int iy, int iz) = 0;
	bool isGeometryLoaded() const;
	
	// return signed distance to closest point on surface (neg. if inside the body).
	virtual double distance(double fx, double fy, double fz) = 0;
	
	void getBounds(double& x1, double& x2, double& y1, double& y2, double& z1, double& z2) const;
	void getCentroid(double p[3]) const;
	void getCenterOfMass(double p[3]) const;
	void getUT(double ut[3]) const;
	void getAngVel(double angVel[3]) const;
	void getMomentOfInertia(double J[6]) const;
	double getDensity() const;
	double getMass() const;
	
	void setScale(double scale);
	void setTranslation(double tX, double tY, double tZ);
	void setBodyIntegrals(double ut[3], double dthetadt[3], double J[6], double mass, double a, double dt);
	void setConstIntegrals(double J[6], double volume);
	
	void moveBody(double dt);
};

#endif /* defined(__CubismUP_3D__GeometryReader__) */
