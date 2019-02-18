//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorPressureGrad_h
#define CubismUP_3D_CoordinatorPressureGrad_h

#include "../SimulationData.h"
#include "GenericCoordinator.h"
#include "GenericOperator.h"

template<int DIRECTION>
class OperatorPressureGradient : public GenericOperator
{
  const double gradPdT;
 public:
  OperatorPressureGradient(double _gradPxdT) : gradPdT(_gradPxdT) { }
  void operator()(const BlockInfo& info, FluidBlock& b) const override
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        if (DIRECTION == 0) b(ix,iy,iz).u += gradPdT;
        if (DIRECTION == 1) b(ix,iy,iz).v += gradPdT;
        if (DIRECTION == 2) b(ix,iy,iz).w += gradPdT;
    }
  }
};

class CoordinatorPressureGradient : public GenericCoordinator
{
 public:
  CoordinatorPressureGradient(SimulationData & s) : GenericCoordinator(s) {}

  void operator()(const double dt)
  {
    check("dp_dx - start");
    const int dir = sim.BCy_flag==wall ? 1 : 2;
    const Real H = sim.extent[dir];
    const Real gradPdt = 8*sim.uMax_forced*sim.nu/H/H * dt;

    const OperatorPressureGradient<0> kernel( gradPdt );
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
        kernel(vInfo[i], *(FluidBlock*)vInfo[i].ptrBlock);
    check("dp_dx - end");
  }

  std::string getName() { return "PressureGradient"; }
};

#endif
