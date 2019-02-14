//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorIC_h
#define CubismUP_3D_CoordinatorIC_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"

class OperatorIC : public GenericOperator
{
 public:
  OperatorIC(const Real u) {}
  ~OperatorIC() {}
  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) block(ix,iy,iz).clear();
  }
};

class OperatorIC_RT : public GenericOperator
{
 public:
  OperatorIC_RT(const Real rhoS) {}
  ~OperatorIC_RT() {}

  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    const Real d = .25, a = .05;//.25
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      Real p[3];
      info.pos(p, ix, iy, iz);
      Real x = p[0] - d*.5;
      Real y = p[1] - d*2.5;
      Real z = p[2] - d*.5;
      //Real eta = -.1*a*cos(2*M_PI*x/d);
      Real eta = -a*.25*std::cos(2*M_PI*x/d)*std::cos(2*M_PI*z/d);
      block(ix,iy,iz).chi = .5+.5*std::tanh((y-eta)/info.h_gridpoint);
    }
  }
};

class OperatorIC_taylorGreen : public GenericOperator
{
  const Real ext[3], uMax;
  const Real a = 2*M_PI / ext[0], b = 2*M_PI / ext[1], c = 2*M_PI / ext[2];
  const Real A = uMax, B = - uMax * ext[1] / ext[0];
 public:
  ~OperatorIC_taylorGreen() {}
  OperatorIC_taylorGreen(const Real extent[3], const Real U): ext{extent[0],extent[1],extent[2]}, uMax(U) {}
  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u = A*std::cos(a*p[0])*std::sin(b*p[1])*std::sin(c*p[2]);
      block(ix,iy,iz).v = B*std::sin(a*p[0])*std::cos(b*p[1])*std::sin(c*p[2]);
    }
  }
};

class OperatorIC_channel : public GenericOperator
{
  const int dir;
  const Real ext[3], uMax, H = ext[dir], FAC = 4*uMax/H/H; // FAC = 0.5*G/mu
  //umax =  0.5*G/mu * 0.25*H*H
 public:
  OperatorIC_channel(const Real extent[3], const Real U, const int _dir):
    dir(_dir), ext{extent[0],extent[1],extent[2]}, uMax(U) {}
  ~OperatorIC_channel() {}

  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u = FAC * p[dir] * (H-p[dir]);
    }
  }
};

class CoordinatorIC : public GenericCoordinator
{
 public:
  CoordinatorIC(SimulationData & s) : GenericCoordinator(s) { }

  template<typename K>
  inline void run(const K kernel) {
    #pragma omp parallel for schedule(static)
    for (size_t i=0; i<vInfo.size(); i++)
      kernel(vInfo[i], *(FluidBlock*)vInfo[i].ptrBlock);
  }
  void operator()(const double dt)
  {
    if(sim.initCond == "zero")
      run(OperatorIC(0));
    if(sim.initCond == "taylorGreen")
      run(OperatorIC_taylorGreen(sim.extent, sim.uMax_forced));
    if(sim.initCond == "channel")
    {
      if(sim.BCx_flag==2) {
        printf("ERROR: channel flow must be periodic or dirichlet in x.\n");
        abort();
      }
      const bool channelY = sim.BCy_flag==2, channelZ = sim.BCz_flag==2;
      if( (channelY && channelZ) or (!channelY && !channelZ) ) {
        printf("ERROR: wrong channel flow BC in y or z.\n");
        abort();
      }
      const int dir = channelY? 1 : 2;
      run(OperatorIC_channel(sim.extent, sim.uMax_forced, dir));
    }
    check("IC - end");
  }

  std::string getName()
  {
    return "IC";
  }
};

#endif
