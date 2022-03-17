//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once
#include "../SimulationData.h"
#include <Cubism/FluxCorrectionMPI.h>

CubismUP_3D_NAMESPACE_BEGIN

class Operator
{
  protected:
   SimulationData & sim;
  public:
   Operator(SimulationData & s) : sim(s) {  }
   virtual ~Operator() = default;
   virtual void operator()(const Real dt) = 0;
   virtual std::string getName() = 0;
};
CubismUP_3D_NAMESPACE_END
