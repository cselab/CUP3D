//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralAnalysis_h
#define CubismUP_3D_SpectralAnalysis_h

#include <vector>
#include <cassert>
#include <cstring>

#include "operators/SpectralManip.h"

CubismUP_3D_NAMESPACE_BEGIN

class SpectralAnalysis
{
 typedef typename FluidGridMPI::BlockType BlockType;
 public:
  // Parameters of the histogram
  int nBin;
  int nyquist;

  // Output of the analysis
  Real tke = 0., eps = 0., tau_integral = 0.;
  Real lambda = 0., uprime = 0., Re_lambda = 0.;
  Real u_avg[3] = {0.};


  Real * k_msr;
  Real * E_msr;

  SpectralAnalysis(SimulationData & s);
  ~SpectralAnalysis();

  void run();
  void dump2File(const int nFile) const;
  void reset();

private:
  SpectralManip * sM;
  void _cub2fftw();
  void _compute();
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralAnalysis_h
