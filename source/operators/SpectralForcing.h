//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Larousslihe.
//

#ifndef CubismUP_3D_SpectralForcing_h
#define CubismUP_3D_SpectralForcing_h

#include "SimulationData.h"
#include "operators/Operator.h"
#include "Cubism/BlockInfo.h"

CubismUP_3D_NAMESPACE_BEGIN

class SpectralManip;

class SpectralForcing : public Operator
{
  Real totalKinEn = 0.0;
  Real viscousDissip = 0.0;
  Real totalKinEnPrev = 0.0;
  Real largeModesKinEn = 0.0;
  void _cub2fftw() const;
  void _fftw2cub(const Real factor) const;

 public:
  SpectralForcing(SimulationData & s);

  void operator()(const double dt);

  std::string getName() { return "SpectralForcing"; }
};

CubismUP_3D_NAMESPACE_END
#endif
