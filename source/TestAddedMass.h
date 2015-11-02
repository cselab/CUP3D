//
//  TestAddedMass.h
//  CubismUP_3D
//
//  Created by Christian Conti on 6/18/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestAddedMass__
#define __CubismUP_3D__TestAddedMass__

#include "Test.h"
#include "GenericCoordinator.h"
#include "Shape.h"

class TestAddedMass : public Test
{
protected:
	ArgumentParser parser;
	
	double nu;
	double minRho;
	double rhoS;
	const int bpd;
	
	Real uBody[3], omegaBody;
	Real gravity[3];
	
	// penalization parameter
	Real lambda;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
	vector<GenericCoordinator *> pipeline;
	
	FluidGrid * grid;
	
	// body
	Shape * shape;
	
	// MPI stuff - required for Hypre
	int rank, nprocs;
	bool bSplit;
	
	// simulation status
	int step, nsteps;
	
	void _ic();
	
public:
	TestAddedMass(const int argc, const char ** argv, const int bpd);
	~TestAddedMass();
	
	void run();
	void check();
};


#endif /* defined(__CubismUP_3D__TestAddedMass__) */
