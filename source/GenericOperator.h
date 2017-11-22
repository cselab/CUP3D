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

class GenericOperator
{
public:
  virtual void operator()(const BlockInfo& info, FluidBlock& block) const = 0;
};

class GenericLabOperator
{
public:
  int stencil_start[3];
  int stencil_end[3];
  StencilInfo stencil;
  // cannot put the templated operator here!
};

#endif
