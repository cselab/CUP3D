//
//  GenericOperator.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/6/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
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
