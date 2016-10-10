//
//  TestMPI.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/18/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef TestMPI_h
#define TestMPI_h

#include <stdio.h>
#include "Test.h"

class TestMPI : public Test
{
private:
	int offset;
	
	string path2file;
	
	void _ic();
	
public:
	TestMPI(const int argc, const char ** argv, const int bpd);
	~TestMPI();
	
	void run();
	void check();
};

#endif /* TestMPI_h */
