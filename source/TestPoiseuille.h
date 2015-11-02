//
//  TestPoiseuille.h
//  CubismUP_3D
//
//  Created by Christian Conti on 5/13/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestPoiseuille__
#define __CubismUP_3D__TestPoiseuille__

#include "Test.h"
#include "GenericCoordinator.h"

class TestPoiseuille : public Test
{
protected:
	double nu;
	double dtCFL, dtLCFL, dtFourier;
	double time, endTime;
	const int bpd;
	int step;
	int rank, nprocs;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	vector<GenericCoordinator *> pipeline;
	
	FluidGrid * grid;
	
	void _ic();
	void _analytical(Real x, Real y, Real z, double t, Real &u, Real &v, Real &w, Real &p);
	
public:
	TestPoiseuille(const int argc, const char ** argv, const int bpd);
	~TestPoiseuille();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestPoiseuille__) */
