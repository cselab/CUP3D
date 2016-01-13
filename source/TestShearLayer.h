//
//  TestShearLayer.h
//  CubismUP_3D
//
//  Created by Christian Conti on 4/28/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestShearLayer__
#define __CubismUP_3D__TestShearLayer__

#include "Test.h"
#include "GenericCoordinator.h"

class TestShearLayer : public Test
{
protected:
	double nu;
	double dtCFL, dtLCFL, dtFourier;
	double time, endTime;
	int step;
	const double rhoS;
	
	// MPI stuff - required for Hypre
	int rank, nprocs;
	
	string path2file;
	
	vector<GenericCoordinator *> pipeline;
	
	void _ic();
	
public:
	TestShearLayer(const int argc, const char ** argv, const int bpd);
	~TestShearLayer();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestShearLayer__) */
