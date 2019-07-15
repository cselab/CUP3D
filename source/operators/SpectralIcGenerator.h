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

#include "operators/SpectralManip.h"

CubismUP_3D_NAMESPACE_BEGIN

class SpectralIcGenerator
{
public:
  typedef typename FluidGridMPI::BlockType BlockType;

  SpectralIcGenerator(SimulationData &s);
  ~SpectralIcGenerator() {}

  void run();

private:
  SpectralManip * sM;

  energySpectrum _generateTarget();
  void _compute();
  void _fftw2cub() const;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralIcGenerator_h
