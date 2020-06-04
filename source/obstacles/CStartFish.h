//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef CubismUP_3D_CStartFish_h
#define CubismUP_3D_CStartFish_h

#include "Fish.h"

CubismUP_3D_NAMESPACE_BEGIN

class CStartFish: public Fish
{
protected:
  Real origC[2] = {(Real)0, (Real)0};
  Real origAng = 0;
public:
  CStartFish(SimulationData&s, cubism::ArgumentParser&p);
  void save(std::string filename = std::string()) override;
  void restart(std::string filename) override;
  void create() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_CStartFish_h
