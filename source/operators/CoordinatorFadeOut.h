//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_CoordinatorFadeOut_h
#define CubismUP_3D_CoordinatorFadeOut_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"

class OperatorFadeOut : public GenericOperator
{
private:
  const Real ext[3], fadeLen[3], iFade[3];
  static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
  inline bool _is_touching(const BlockInfo& i) const
  {
    Real maxP[3], minP[3]; i.pos(minP, 0, 0, 0);
    i.pos(maxP, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1);
    const bool touchW= fadeLen[0]>=minP[0], touchE= fadeLen[0]>=ext[0]-maxP[0];
    const bool touchS= fadeLen[1]>=minP[1], touchN= fadeLen[1]>=ext[1]-maxP[1];
    const bool touchB= fadeLen[2]>=minP[2], touchF= fadeLen[2]>=ext[2]-maxP[2];
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }
  inline Real fade(const BlockInfo&i, const int x,const int y,const int z) const
  {
    Real p[3]; i.pos(p, x, y, z);
    const Real zt = iFade[2] * std::max(Real(0), fadeLen[2] -(ext[2]-p[2]) );
    const Real zb = iFade[2] * std::max(Real(0), fadeLen[2] - p[2] );
    const Real yt = iFade[1] * std::max(Real(0), fadeLen[1] -(ext[1]-p[1]) );
    const Real yb = iFade[1] * std::max(Real(0), fadeLen[1] - p[1] );
    const Real xt = iFade[0] * std::max(Real(0), fadeLen[0] -(ext[0]-p[0]) );
    const Real xb = iFade[0] * std::max(Real(0), fadeLen[0] - p[0] );
    return 1-std::pow(std::min( std::max({zt,zb,yt,yb,xt,xb}), (Real)1), 2);
  }

public:
  OperatorFadeOut(const Real buf[3], const Real extent[3]) :
  ext{extent[0],extent[1],extent[2]}, fadeLen{buf[0],buf[1],buf[2]},
  iFade{1/(buf[0]+EPS), 1/(buf[1]+EPS), 1/(buf[2]+EPS)} {}

  void operator()(const BlockInfo& info, FluidBlock& b) const
  {
    if( _is_touching(info) )
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const Real FADE = fade(info, ix, iy, iz);
      b(ix,iy,iz).u *= FADE; b(ix,iy,iz).v *= FADE; b(ix,iy,iz).w *= FADE;
    }
  }
};

class CoordinatorFadeOut : public GenericCoordinator
{
 public:
    CoordinatorFadeOut(SimulationData & s) : GenericCoordinator(s) { }

  void operator()(const double dt)
  {
    check((std::string)"FadeOut - start");

    const int N = vInfo.size();
    #pragma omp parallel
    {
      OperatorFadeOut kernel(sim.fadeOutLengthU, sim.extent);
      #pragma omp for schedule(static)
      for (int i=0; i<N; i++) {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        kernel(info, b);
      }
    }
    check((std::string)"FadeOut - end");
  }

  std::string getName()
  {
    return "FadeOut";
  }
};

#endif
