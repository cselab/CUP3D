//
//  Test.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_Test_h
#define CubismUP_3D_Test_h

#include <stdio.h>
#include "Definitions.h"
#include "CoordinatorVorticity.h"

class Test
{
protected:
	ArgumentParser parser;
	
	const int bpdx, bpdy, bpdz;
	int rank, nprocs;
	int nprocsx, nprocsy, nprocsz;
	
	FluidGridMPI * grid;
	
	vector<BlockInfo> vInfo;
	
public:
	Test(const int argc, const char ** argv, const int bpd) : parser(argc,argv), bpd(bpd)
	{
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

		bpdx = parser("-bpdx").asInt();
		bpdy = parser("-bpdy").asInt();
		bpdz = parser("-bpdz").asInt();
		nprocsx = parser("-nprocsx").asInt(nprocs);
		nprocsy = parser("-nprocsy").asInt(1);
		nprocsz = parser("-nprocsz").asInt(1);
		assert(bpdx%nprocsx == 0 && bpdy%nprocsy == 0 && bpdz%nprocsz == 0);
		cout << "nranks: " << nprocsx << " " << nprocsy << " " << nprocsz << endl;
		cout << "bpd/rank: " << bpdx/nprocsx << " " << bpdy/nprocsy << " " << bpdz/nprocsz << endl;
		grid = new FluidGridMPI(nprocsx, nprocsy, nprocsz, bpdx/nprocsx, bpdy/nprocsy, bpdz/nprocsz);
		vInfo = grid->getBlocksInfo();
	}
	
	virtual ~Test()
	{
		delete grid;
	}
	
	virtual void run() = 0;
	virtual void check() = 0;
};

#endif
