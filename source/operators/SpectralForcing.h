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

class SpectralForcing : public Operator
{
 public:
  SpectralForcing(SimulationData & s) : Operator(s) {}

  void operator()(const double dt);

  std::string getName() { return "SpectralForcing"; }
};

CubismUP_3D_NAMESPACE_END
#endif
