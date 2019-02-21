//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "accfft_utils.h"
#include "accfft_common.h"

#ifndef CUP_SINGLE_PRECISION
  #include "accfft_gpu.h"
  typedef double acc_c[2];
  typedef accfft_plan_gpu acc_plan;
  #define accfft_local_size accfft_local_size_dft_r2c_gpu
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpu
  #define accfft_delplan accfft_destroy_plan_gpu
  #define accfft_clean accfft_cleanup_gpu
  #define accfft_exec_r2c accfft_execute_r2c_gpu
  #define accfft_exec_c2r accfft_execute_c2r_gpu
#else
  #include "accfft_gpuf.h"
  typedef float acc_c[2];
  typedef accfft_plan_gpuf acc_plan;
  #define accfft_local_size accfft_local_size_dft_r2c_gpuf
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpuf
  #define accfft_delplan accfft_destroy_plan_gpu
  #define accfft_clean accfft_cleanup_gpuf
  #define accfft_exec_r2c accfft_execute_r2c_gpuf
  #define accfft_exec_c2r accfft_execute_c2r_gpuf
#endif
