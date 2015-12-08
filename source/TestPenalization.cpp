//
//  TestPenalization.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 2/3/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "TestPenalization.h"
#include <sstream>
#include <cmath>

#include "CoordinatorIC.h"
#include "CoordinatorPenalization.h"

void TestPenalization::_ic()
{
	Real center[3] = {.5,.5,.5};
	Real radius = .1;
	
	shape = new Sphere(center, radius, (Real).1, (Real)2, (Real)2);
	
	CoordinatorIC coordIC(shape,1.,grid);
	coordIC(0);
}

TestPenalization::TestPenalization(const int argc, const char ** argv, const int bpd, const double dt) : Test(argc, argv), bpd(bpd), lambda(1e6), uBody{0,0.5,0.1}, dt(dt)
{
	grid = new FluidGrid(bpd,bpd,bpd);
	
	// setup initial condition
	
	// output settings
	path2file = parser("-file").asString("../data/testPenalization");
	_ic();
}

TestPenalization::~TestPenalization()
{
	delete grid;
	delete shape;
}

void TestPenalization::run()
{
	vector<BlockInfo> vInfo = grid->getBlocksInfo();
	
	Real omegaBody = 0;
	CoordinatorPenalization coordPenalization(&uBody[0], &uBody[1], &uBody[2], shape, &lambda, grid);
	
	for (int i=0; i<10; i++)
	{
		coordPenalization(dt);
	
		stringstream ss;
		ss << path2file << "-bpd" << bpd << "-step" << i << ".vti";
	
		dumper.Write(*grid, ss.str());
	}
}

void TestPenalization::check()
{
}
