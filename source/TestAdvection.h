//
//  TestAdvection.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/8/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestAdvection__
#define __CubismUP_3D__TestAdvection__

#include <stdio.h>
#include "Test.h"
#include "Layer.h"

class TestAdvection : public Test
{
private:
	double time, dt;
	int bpd;
	int testCase; // 0: linear, 1: rigid rotation
	int nsteps;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	FluidGrid * grid;
	
	Layer * vorticityIC;
	
	void _icLinear();
	void _icVortex();
	void _icBurger();
	
public:
	TestAdvection(const int argc, const char ** argv, int testCase, const int bpd, const double dt, const int nsteps);
	~TestAdvection();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestAdvection__) */
