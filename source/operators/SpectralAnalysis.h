//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralAnalysis_h
#define CubismUP_3D_SpectralAnalysis_h

#include "../SimulationData.h"

#include <Cubism/BlockInfo.h>

#include <cassert>
#include <cstring>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

class SpectralManip;

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
  Real u_avg[3] = {(Real)0, (Real)0, (Real)0};
  Real unorm_avg = 0;
  Real * k_msr;
  Real * E_msr;

  bool bComputeCs2Spectrum = false;
  Real * cs2_msr;

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
