//
//  TestGeometry.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/4/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef TestGeometry_h
#define TestGeometry_h

#include <stdio.h>
#include "Test.h"
#include "Shape.h"

class TestGeometry : public Test
{
private:
	Shape * shape;
	const int testCase; // 0: forced, 1: from flow
	double dt;
	
	string path2file;
	
	void _ic();
	
public:
	TestGeometry(const int argc, const char ** argv, const int bpd);
	~TestGeometry();
	
	void run();
	void check();
};

#endif /* TestGeometry_h */
