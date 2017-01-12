//
//  Simulation_Fluid.h
//  CubismUP_2D
//
//	Base class for fluid simulations from which any fluid simulation case should inherit
//	Contains the base structure and interface that any fluid simulation class should have
//
//  Created by Christian Conti on 3/25/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_2D_Simulation_Fluid_h
#define CubismUP_2D_Simulation_Fluid_h

//#include "Definitions.h"
#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "ProcessOperatorsOMP.h"

#include "CoordinatorIC.h"
#include "CoordinatorVorticity.h"
#include "CoordinatorAdvection.h"
#include "CoordinatorDiffusion.h"
#include "CoordinatorPenalization.h"
#include "CoordinatorComputeShape.h"
#include "CoordinatorPressure.h"
#include "CoordinatorFadeOut.h"
#include "IF3D_ObstacleVector.h"
#include "IF3D_ObstacleFactory.h"
//#include "IF3D_StefanFishOperator.h"
//#include "IF3D_CarlingFishOperator.h"
//#include "IF3D_SphereObstacleOperator.h"
//#include "IF3D_ForcedSphereObstacleOperator.h"
#ifdef _USE_LZ4_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif

#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace ComputationDiagnostics
{
    static void print_memory_usage(long & peak_rss_bytes, long & current_rss_bytes)
    {
        /*
        peak_rss_bytes = -1;
        current_rss_bytes=-1;

        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);

        peak_rss_bytes = rusage.ru_maxrss*1024;
        //printf("peak resident set size = %ld bytes (%.2lf Mbytes)\n", peak_rss_bytes, peak_rss_bytes/(1024.0*1024.0));

        long rss = 0;
        FILE* fp = NULL;
        if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL ) {
            return;
        }

        if ( fscanf( fp, "%*s%ld", &rss ) != 1 ) {
            fclose( fp );
            return;
        }
        fclose( fp );

        current_rss_bytes = rss * sysconf( _SC_PAGESIZE);
        //printf("current resident set size = %ld bytes (%.2lf Mbytes)\n", current_rss_bytes, current_rss_bytes/(1024.0*1024.0));
        */
    }
}

class Simulation
{
protected:
	ArgumentParser parser;
	Profiler profiler;
	Communicator * communicator;
#ifdef _USE_LZ4_
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, FluidVPStreamer> waveletdumper_grid;
#endif

	// grid
	int rank, nprocs;
	int nprocsx, nprocsy, nprocsz;
	int bpdx, bpdy, bpdz;

	// simulation status
	int step, nsteps;
	Real dt, time, endTime, dtCFL, dtFourier;

	// simulation settings
    Real uinf[3], re, nu, length, CFL, lambda;
    bool bDump, bRestart, bDLM, verbose, b2Ddump;

	// output
	int dumpFreq, saveFreq;
	Real dumpTime, saveTime, saveClockPeriod, maxClockDuration;
	string path2file, path4serialization;

	FluidGridMPI * grid;
    vector<BlockInfo> vInfo;
	//The protagonist
    IF3D_ObstacleVector* obstacle_vector;
    //The antagonist
	vector<GenericCoordinator*> pipeline;

    void areWeDumping(Real & nextDumpTime);
    void _serialize(Real & nextSaveTime);
    void _dump(const string append);
    void _deserialize();

    void parseArguments();
    void setupObstacles();
    void setupOperators();
    void _selectDT();
    void setupGrid();
    void _ic();

public:
    Simulation(const int argc, char ** argv, Communicator* comm = nullptr) :
    parser(argc,argv), communicator(comm), rank(0), nprocs(1), nprocsx(-1), nprocsy(-1),
	nprocsz(-1), bpdx(-1), bpdy(-1), bpdz(-1), step(0), nsteps(0), dt(0), time(0), endTime(0),
	dtCFL(0), dtFourier(0), uinf{0.0, 0.0, 0.0}, re(0), nu(0), length(0), CFL(0), lambda(0),
	bDump(false), bRestart(false), bDLM(false), verbose(false), b2Ddump(false),
	dumpFreq(0), saveFreq(0), dumpTime(0), saveTime(0), saveClockPeriod(0), maxClockDuration(1e9)
	{
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		char hostname[1024];
		hostname[1023] = '\0';
		gethostname(hostname, 1023);
		const int nthreads = omp_get_max_threads();
		printf("Rank %d (of %d) with %d threads on host Hostname: %s\n", rank, nprocs, nthreads, hostname);
	}

	virtual ~Simulation()
	{
		delete grid;
		while(!pipeline.empty()) {
			GenericCoordinator * g = pipeline.back();
			pipeline.pop_back();
			delete g;
		}
	}

    virtual void init();

    virtual void simulate();
};

#endif
