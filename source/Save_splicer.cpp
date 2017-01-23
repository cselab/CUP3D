//
//  Save_splicer_Fluid.h
//  CubismUP_2D
//
//	Base class for fluid Save_splicers from which any fluid Save_splicer case should inherit
//	Contains the base structure and interface that any fluid Save_splicer class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "Save_splicer.h"
#include <HDF5Dumper_MPI.h>
#include "ProcessOperatorsOMP.h"
#include <chrono>

void Save_splicer::setupGrid()
{
	parser.set_strict_mode();
	bpdx = parser("-bpdx").asInt();
	bpdy = parser("-bpdy").asInt();
	bpdz = parser("-bpdz").asInt();
	nprocsx = parser("-nprocsx").asInt();
	parser.unset_strict_mode();
	nprocsy = parser("-nprocsy").asInt(1);
	nprocsz = parser("-nprocsz").asInt(1);
	assert(bpdx%nprocsx==0 && bpdy%nprocsy==0 && bpdz%nprocsz==0);
	if(rank==0) {
		printf("Blocks per dimension: [%d %d %d]\n",bpdx,bpdy,bpdz);
		printf("Nranks per dimension: [%d %d %d]\n",nprocsx,nprocsy,nprocsz);
		printf("Resident blocks per dim: [%d %d %d]\n",bpdx/nprocsx,bpdy/nprocsy,bpdz/nprocsz);
	}
	bpdx /= nprocsx;
	bpdy /= nprocsy;
	bpdz /= nprocsz;
	grid = new FluidGridMPI(nprocsx, nprocsy, nprocsz, bpdx, bpdy, bpdz);
	assert(grid != NULL);
    vInfo = grid->getBlocksInfo();
}

void Save_splicer::parseArguments()
{
    path2list = parser("-list").asString("list.dat");
    path4serialization = parser("-serialization").asString("./");
    path4deserialization  = parser("-deserialization").asString("./");
}

void Save_splicer::load_and_dump(string path2file)
{
	string restartfile = path4deserialization+"/"+path2file;
	FILE * f = fopen(restartfile.c_str(), "r");
	if (f == NULL) {
		printf("Could not restart....\n");
		return;
	}
	assert(f != NULL);

	float val = -1;
	fscanf(f, "time: %e\n", &val);
	assert(val>=0);
	time=val;

	int step_id_fake = -1;
	fscanf(f, "stepid: %d\n", &step_id_fake);
	assert(step_id_fake >= 0);
	step = step_id_fake;

	int ret = 0;
	ret = fscanf(f, "uinfx: %e\n", &val);
	if (ret) uinf[0] = val;
	ret = fscanf(f, "uinfy: %e\n", &val);
	if (ret) uinf[1] = val;
	ret = fscanf(f, "uinfz: %e\n", &val);
	if (ret) uinf[2] = val;

	fclose(f);
	printf("DESERIALIZATION: time is %f and step id is %d\n", time, (int)step);

	MPI_Barrier(MPI_COMM_WORLD);
	{
	stringstream ssR;
	ssR<<"restart_"<<std::setfill('0')<<std::setw(9)<<step;
	if(rank==0) cout<<"Restart from "<<path4deserialization<<" "<<ssR.str()<<endl;
	fflush(0);
	#if defined(_USE_HDF_)
	ReadHDF5_MPI<FluidGridMPI,StreamerHDF5>(*grid,ssR.str(),path4deserialization);
	#endif
	}
	MPI_Barrier(MPI_COMM_WORLD);
	{
	stringstream ssR;
	ssR<<"restart_"<< std::setfill('0')<<std::setw(9)<<step;
	if (rank==0)  cout << "Saving into "<<path4serialization<<" "<< ssR.str() << endl;
	fflush(0);
	#if defined(_USE_HDF_)
	//DumpHDF5_MPI_Channel<StreamerHDF5,0>(*grid, time, ssR.str(), path4serialization);
	#endif
	}
}

void Save_splicer::run()
{
    parseArguments();
    setupGrid();

		ifstream in(path2list.c_str());
		std::string path2file, line;
		if(in.good())
		while (getline(in, line)) {
			istringstream line_in(line);
			while(line_in >> path2file)
				load_and_dump(path2file);
		}
}
