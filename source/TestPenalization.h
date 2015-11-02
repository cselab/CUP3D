//
//  TestPenalization.h
//  CubismUP_3D
//
//  Created by Christian Conti on 2/3/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestPenalization__
#define __CubismUP_3D__TestPenalization__

#include <stdio.h>
#include "Test.h"
#include "Shape.h"

class TestPenalization : public Test
{
private:
	int bpd;
	double dt;
	Real lambda;
	Real uBody[3];
	Shape * shape;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	FluidGrid * grid;
	
	void _ic();
	
public:
	TestPenalization(const int argc, const char ** argv, const int bpd, const double dt);
	~TestPenalization();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestPenalization__) */
