//
//  CTReader.h
//  CubismUP_3D
//
//	Reads a series of CT slices and computes the signed distance function
//
//	Created by Christian Conti on 08/15/13. (MPCF)
//  Created by Christian Conti on 11/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__CTReader__
#define __CubismUP_3D__CTReader__

#pragma once

#include <fstream>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <cassert>
#include <iomanip>

#include "GeometryHelpers.h"
#include "GeometryReader.h"

class GeometryReaderCT : public GeometryReader
{
protected:
	int subsampling;
	int sx, sy, sz;
	int buffer;
	bool binary;
	int nchannels;
	
	virtual void readFromFile(const int slice, const int startslice);
	virtual void parse(string filename);
	
public:
	// this is the grid on which to interpolate the mesh,
	// it could be eventually placed outside of this class and
	// be a cubism grid
	Geometry::Grid density;
	
	// constructor loading OBJ file
	GeometryReaderCT(const string filename, Geometry::Properties & properties, int gridsize, const Real scaleFactor, const Geometry::Point transFactor);
	
	// clean up
	~GeometryReaderCT();
	
	virtual void load(const string filename);
	
	virtual void sdf(const Real scaleFactor=1, const Geometry::Point transFactor=Geometry::Point(), double comx=1, double comy=1, double comz=1);
	
	// check if a point lies inside the body
	virtual int isPointInside(int ix, int iy, int iz);
	
	// return signed distance to closest point on surface (neg. if inside the body).
	virtual double distance(double fx, double fy, double fz);
};


#endif /* defined(__CubismUP_3D__CTReader__) */
