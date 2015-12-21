//
//  WavefrontReader.h
//  CubismUP_3D
//
//  Reads a wavefront object file and computes the signed distance function
//
//	Created by Christian Conti on 06/14/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__WavefrontReader__
#define __CubismUP_3D__WavefrontReader__

#pragma once

#include <fstream>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <cassert>

#include <HDF5Dumper.h>

#include "GeometryReader.h"

typedef ScalarGrid GScalar;

class GeometryReaderOBJ : public GeometryReader
{
protected:
	vector<Geometry::Point> vertices;
	vector<Geometry::Point> vnormals;
	vector<Geometry::Triangle> triangles;
	
	double grid_maxsize;
	GScalar * cgrid;
	// grid lut
	//	0: outside
	//	1: boundary
	//	2: inside
	//	3: unset
	int * gridlut;
	
	const Real isosurface;
	
	void parse(string filename);
	
public:
	// constructor loading OBJ file	
	GeometryReaderOBJ(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor, const Real isosurface);
	
	// clean up
	~GeometryReaderOBJ();
	
	virtual void load(const string filename);
	
	virtual void sdf();
	
	// check if a point lies inside the body
	virtual int isPointInside(int ix, int iy, int iz);
	
	// return signed distance to closest point on surface (neg. if inside the body).
	virtual double distance(double fx, double fy, double fz);
};


#endif /* defined(__CubismUP_3D__WavefrontReader__) */
