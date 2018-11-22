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
  const Real extent[3];
  const Real buffer;

  inline bool _is_touching(const BlockInfo& info, const Real h) const
  {
    Real max_pos[3],min_pos[3];
    info.pos(min_pos, 0, 0, 0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    #ifdef BC_PERIODICZ
      const bool touchF = false;
      const bool touchB = false;
    #else
      const bool touchF = h + buffer >= extent[2] - max_pos[2];
      const bool touchB = h + buffer >= min_pos[2];
    #endif
    const bool touchE = h + buffer >= extent[0] - max_pos[0];
    const bool touchW = h + buffer >= min_pos[0];

    const bool touchS = h + buffer >= min_pos[1];
    const bool touchN = h + buffer >= extent[1] - max_pos[1];
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }

public:
  OperatorFadeOut(const Real buf, const Real ext[3]): extent{ext[0],ext[1],ext[2]},buffer(buf){}

  void operator()(const BlockInfo& info, FluidBlock& b) const
  {
    const Real h = info.h_gridpoint, iWidth = 1/buffer;
    if(_is_touching(info,h))
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      Real p[3];
      info.pos(p, ix, iy, iz);
      #ifdef BC_PERIODICZ
        const Real dzt = 0, dzb = 0;
      #else
        const Real dzt = std::max(Real(0), h+buffer -(extent[2]-p[2]) );
        const Real dzb = std::max(Real(0), h+buffer -p[2] );
      #endif
      const Real dxt = std::max(Real(0), h+buffer -(extent[0]-p[0]) );
      const Real dyt = std::max(Real(0), h+buffer -(extent[1]-p[1]) );
      const Real dxb = std::max(Real(0), h+buffer -p[0] );
      const Real dyb = std::max(Real(0), h+buffer -p[1] );
      // max distance in killing zone 0 <= out <= h+buffer
      const Real out = max(max(max(dxt,dxb), max(dyt,dyb)), max(dzt,dzb));
      // 1 at buffer start, 0 at buffer end (2 grid points before border)
      const Real fade = std::max(Real(0), std::cos(Real(.5*M_PI*out*iWidth)));
      // smooth within killing zone (factor <= 1) and kill at very boundaries (factor < 0)
      b(ix,iy,iz).u = b(ix,iy,iz).u*fade;
      b(ix,iy,iz).v = b(ix,iy,iz).v*fade;
      b(ix,iy,iz).w = b(ix,iy,iz).w*fade;
      //b(ix,iy,iz).p = b(ix,iy,iz).p*fade;
    }
  }
};

class CoordinatorFadeOut : public GenericCoordinator
{
protected:
  const Real buffer;
public:
    CoordinatorFadeOut(FluidGridMPI * g, const Real _buffer = 0)
  : GenericCoordinator(g), buffer(_buffer)
  { }

  void operator()(const double dt)
  {
    check((std::string)"FadeOut - start");

    const int N = vInfo.size();
    const Real h = grid->getH();
    const Real ext[3] = {
        h*grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        h*grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        h*grid->getBlocksPerDimension(2)*FluidBlock::sizeZ
    };
    const Real B = std::max(3 * h, buffer);
    #pragma omp parallel
    {
      OperatorFadeOut kernel(B, ext);
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
