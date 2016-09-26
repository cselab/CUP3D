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

/*
namespace ComputationDiagnostics
{
    static void print_memory_usage(long & peak_rss_bytes, long & current_rss_bytes)
    {
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
    }
}

namespace FlowDiagnostics
{
    struct diagContainer
    {
        // three linear invariants (conserved in inviscid and viscous flows)
        Real circ; // conservation of vorticity: int w dx = 0
        Real linImpulse[2]; // conservation of linear impulse: int u dx = int (x cross w) dx
        Real angImpulse; // conservation of angular impulse: int (x cross u)dx = 1/3 int (x cross (x cross w) )dx
        
        // two more important: maxvor and enstrophy
        Real ens,maxvor;
        
        diagContainer()
        {*this=0;}
        
        diagContainer & operator=(const Real & rhs)
        {
            circ = rhs;
            linImpulse[0] = rhs;
            linImpulse[1] = rhs;
            angImpulse = rhs;
            ens=rhs;
            maxvor=rhs;
            return *this;
        }
        
        diagContainer & operator+=(const diagContainer & rhs)
        {
            circ += rhs.circ;
            linImpulse[0] += rhs.linImpulse[0];
            linImpulse[1] += rhs.linImpulse[1];
            angImpulse += rhs.angImpulse;
            ens += rhs.ens;
            maxvor = std::max(maxvor,rhs.maxvor);
            
            return *this;
        }
    };
    
    struct basicDiagnostics
    {
        diagContainer data;
        
        basicDiagnostics()
        {data=Real(0);}
        
        // split constructor (dont need copy constructor)
        basicDiagnostics(const basicDiagnostics& c, tbb::split):vInfo(c.vInfo), coll(c.coll)
        {data=Real(0);}
        
        void join(const basicDiagnostics& j)
        {
            data+=j.data;
        }
        
        
        void operator()(const BlockInfo& info, FluidBlock& block) const
        {
            // this implementation considers that the Euler updates has already happened
            // do we need a finite state machine coordinating operators?
            diagContainer blockData;
            blockData=Real(0);
            This works only for 2D and this text ensures that you will not compile this D:
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    Real x[2] = {0,0};
                    info.pos(x,ix,iy);
                    
                    const Real w = block(ix,iy).omega;
                    const Real u[2] = { block(ix,iy).u,
                                        block(ix,iy).v };
                    
                    blockData.circ += w;
                    blockData.linImpulse[0] += x[1]*w;
                    blockData.linImpulse[1] -= x[0]*w;
                    blockData.angImpulse -= (x[0]*x[0]+x[1]*x[1])*w;
                    
                    blockData.ens += w*w;
                    blockData.maxvor = std::max(blockData.maxvor,w);
                }
            const Real dA = info.h[0]*info.h[1];
            
            blockData.circ*=dA;
            blockData.linImpulse[0]*=dA;
            blockData.linImpulse[1]*=dA;
            blockData.angImpulse*=dA;
            blockData.ens*=dA;
            
            data+=blockData;
        }
    };
}
*/

class Simulation
{
protected:
	ArgumentParser parser;
	Profiler profiler;
#ifdef _USE_LZ4_
	SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, FluidVPStreamer> waveletdumper_grid;
#endif
	// Serialization
	bool bPing; // needed for ping-pong scheme
	string path4serialization;
	
	//The protagonist
    IF3D_ObstacleVector* obstacle_vector;
    //The antagonist
	vector<GenericCoordinator*> pipeline;
	Communicator * communicator;

	// grid
	int rank, nprocs;
	int nprocsx, nprocsy, nprocsz;
	int bpdx, bpdy, bpdz;
	FluidGridMPI * grid;
    vector<BlockInfo> vInfo;
	
	// simulation status
	int step, nsteps;
	double dt, time, endTime;
    double uinf[3], uinf_dummy[3], re, nu, length;
    double dtCFL, dtLCFL, dtFourier;
	
	// simulation settings
	double CFL, LCFL, lambda, theta;
    bool bDump, bRestart, bDLM, verbose, b2Ddump;
	
	// output
	int dumpFreq, saveFreq;
	double dumpTime, saveTime;
    //double nextDumpTime, nextSaveTime;
	string path2file;
	//SerializerIO_ImageVTK<FluidGrid, FluidVTKStreamer> dumper;
	
    void areWeDumping(double & nextDumpTime);
    void _serialize(double & nextSaveTime);
    void _dump(const string append);
    void _deserialize();

    void parseArguments();
    void setupObstacles();
    void setupOperators();
    void _selectDT();
    void setupGrid();
    void _ic();
	
public:
    Simulation(const int argc, char ** argv, Communicator* comm) :
    parser(argc,argv), bpdx(-1), bpdy(-1), step(0), nsteps(0), endTime(0), time(0), dt(0),
	rank(0), nprocs(1), bPing(false), uinf{0.0,0.0}, uinf_dummy{0.0,0.0}, re(0), length(0), nu(0),
	dtCFL(0), dtLCFL(0), dtFourier(0), CFL(0), LCFL(0), lambda(0), bDump(false), bRestart(false), b2Ddump(false),
	verbose(false), bDLM(false), dumpFreq(0), saveFreq(0), dumpTime(0), saveTime(0), communicator(comm)
	{
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
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
