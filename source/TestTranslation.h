//
//  TestTranslation.h
//  CubismUP_3D
//
//	A neutrally buoyant disk immersed in a fluid with u_inf should move at the same velocity
//
//  Created by Christian Conti on 3/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestTranslation__
#define __CubismUP_3D__TestTranslation__

#include <stdio.h>
#include "Test.h"
#include "Shape.h"

class TestTranslation : public Test
{
private:
	int bpd;
	Real uBody[3];
	Shape * shape;
	const int testCase; // 0: forced, 1: from flow
	double dt;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	FluidGrid * grid;
	
	void _ic();
	
public:
	TestTranslation(const int argc, const char ** argv, const int testCase, const int bpd, const double dt);
	~TestTranslation();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestTranslation__) */
