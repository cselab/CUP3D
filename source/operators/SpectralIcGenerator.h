//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralIcGenerator_h
#define CubismUP_3D_SpectralIcGenerator_h

#include <vector>
#include <cassert>
#include <cstring>

#include "../SimulationData.h"
#include "Cubism/BlockInfo.h"

CubismUP_3D_NAMESPACE_BEGIN

class SpectralManip;

class SpectralIcGenerator
{
public:
  typedef typename FluidGridMPI::BlockType BlockType;

  SpectralIcGenerator(SimulationData &s);
  ~SpectralIcGenerator() {}

  void run();

private:
  SpectralManip * sM;

  void _generateTarget(std::vector<Real>& k, std::vector<Real>& E);
  void _fftw2cub() const;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralIcGenerator_h
