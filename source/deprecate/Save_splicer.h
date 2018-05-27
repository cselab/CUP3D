//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//


#ifndef CubismUP_2D_Save_splicer_Fluid_h
#define CubismUP_2D_Save_splicer_Fluid_h

//#include "Definitions.h"
#include "GenericOperator.h"
#include "GenericCoordinator.h"

#ifdef _USE_LZ4_
#include "SerializerIO_WaveletCompression_MPI_Simple.h"
#endif

#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


class Save_splicer
{
protected:
  ArgumentParser parser;
  #ifdef _USE_LZ4_
  //SerializerIO_WaveletCompression_MPI_SimpleBlocking<FluidGridMPI, FluidVPStreamer> waveletdumper_grid;
  #endif

  // grid
  int rank, nprocs;
  int nprocsx, nprocsy, nprocsz;
  int bpdx, bpdy, bpdz;

  // Save_splicer status
  int step;
  Real time;

  // Save_splicer settings
  Real uinf[3];
  string path2list, path4serialization, path4deserialization;

  FluidGridMPI * grid;
  vector<BlockInfo> vInfo;

  void load_and_dump(string path2file);
  void parseArguments();
  void setupGrid();

public:
  Save_splicer(const int argc, char ** argv) :
  parser(argc,argv), rank(0), nprocs(1), nprocsx(-1), nprocsy(-1),
  nprocsz(-1), bpdx(-1), bpdy(-1), bpdz(-1), step(0), time(0), uinf{0,0,0}
  {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    char hostname[1024];
    hostname[1023] = '\0';
    gethostname(hostname, 1023);
    const int nthreads = omp_get_max_threads();
    printf("Rank %d (of %d) with %d threads on host Hostname: %s\n", rank, nprocs, nthreads, hostname);
  }

  virtual ~Save_splicer()
  {
    delete grid;
  }

  virtual void run();
};

#endif
