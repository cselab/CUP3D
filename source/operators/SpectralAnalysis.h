//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralAnalysis_h
#define CubismUP_3D_SpectralAnalysis_h

#include "SimulationData.h"
#include "Cubism/BlockInfo.h"

#include <vector>
#include <cassert>
#include <cstring>

CubismUP_3D_NAMESPACE_BEGIN

class SpectralManip;

class SpectralAnalysis
{
 typedef typename FluidGridMPI::BlockType BlockType;
 public:

  Real u_avg[3] = {(Real)0, (Real)0, (Real)0};
  Real unorm_avg = 0;

  SpectralAnalysis(SimulationData & s);
  ~SpectralAnalysis();

  void run();
  void dump2File(const int nFile) const;
  void reset();

private:
  SpectralManip * sM;
  void _cub2fftw();
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralAnalysis_h
