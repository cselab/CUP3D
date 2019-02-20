//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include <fftw3-mpi.h>
#ifndef CUP_SINGLE_PRECISION
#define _FFTW_(s) fftw_##s
typedef fftw_complex fft_c;
typedef fftw_plan fft_plan;
#else
#define _FFTW_(s) fftwf_##s
typedef fftwf_complex fft_c;
typedef fftwf_plan fft_plan;
#endif /* CUP_SINGLE_PRECISION */
