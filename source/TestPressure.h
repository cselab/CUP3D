//
//  TestPressure.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestPressure__
#define __CubismUP_3D__TestPressure__

#include <stdio.h>
#include "Test.h"

class TestPressure : public Test
{
private:
	int bpd;
	const int ic, solver;
	double dt;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	FluidGrid * grid;
	int rank, nprocs;
	
	void _ic();
	
public:
	TestPressure(const int argc, const char ** argv, const int solver, const int ic, const int bpd, const double dt);
	~TestPressure();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestPressure__) */
