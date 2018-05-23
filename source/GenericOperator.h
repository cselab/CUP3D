//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_Operator_h
#define CubismUP_3D_Operator_h
#include "Definitions.h"

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
