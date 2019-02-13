//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_Operator_h
#define CubismUP_3D_Operator_h

#include "../Definitions.h"

#include <array>

class GenericOperator
{
public:
  virtual void operator()(const BlockInfo& info, FluidBlock& block) const = 0;
};

class GenericLabOperator
{
public:
  std::array<int, 3> stencil_start;
  std::array<int, 3> stencil_end;
  StencilInfo stencil;
  // cannot put the templated operator here!
};

#endif
