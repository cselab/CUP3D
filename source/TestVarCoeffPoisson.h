//
//  TestVarCoeffPoisson.h
//  CubismUP_3D
//
//  Created by Christian Conti on 1/23/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__TestVarCoeffPoisson__
#define __CubismUP_3D__TestVarCoeffPoisson__

#include <stdio.h>
#include "Test.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif // _MULTIGRID_

class TestVarCoeffPoisson : public Test
{
private:
	const int ic;
	const int bpd;
	
	string path2file;
	SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
#ifdef _MULTIGRID_
	MultigridHypre mg;
#endif // _MULTIGRID_
	
	FluidGrid * grid;
	
	void _ic();
	
public:
	TestVarCoeffPoisson(const int argc, const char ** argv, const int ic, const int bpd);
	~TestVarCoeffPoisson();
	
	void run();
	void check();
};

#endif /* defined(__CubismUP_3D__TestVarCoeffPoisson__) */
