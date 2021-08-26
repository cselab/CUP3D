//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_CarlingFish_h
#define CubismUP_3D_CarlingFish_h

#include "Fish.h"

CubismUP_3D_NAMESPACE_BEGIN

class CarlingFishMidlineData;

class CarlingFish: public Fish
{
 public:
  CarlingFish(SimulationData&s, cubism::ArgumentParser&p);
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_CarlingFish_h
