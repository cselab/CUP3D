//
//  common.h
//  CubismUP_3D
//
//  Created by Christian Conti on 11/2/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_common_h
#define CubismUP_3D_common_h

#include <cassert>
#include <sstream>
#include <cmath>
#include <cstdio>

// utmost import to be defined before including cubism

#ifndef _SP_COMP_
typedef double Real;
#else // _SP_COMP_
typedef float Real;
#endif // _SP_COMP_

#include <mpi.h>
#include <omp.h>

//this is all cubism file we need
#include <ArgumentParser.h>
#include <Grid.h>
#include <GridMPI.h>
#include <BlockInfo.h>
#ifdef _VTK_
#include <SerializerIO_ImageVTK.h>
#endif
#include <HDF5Dumper_MPI.h>
#include <ZBinDumper_MPI.h>
#include <BlockLab.h>
#include <BlockLabMPI.h>
#include <Profiler.h>
#include <StencilInfo.h>

#include "Timer.h"

#endif
