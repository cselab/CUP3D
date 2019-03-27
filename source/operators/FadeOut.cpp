//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/FadeOut.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

class KernelFadeOut
{
 private:
  const Real ext[3], fadeLen[3], iFade[3];
  static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
  inline bool _is_touching(const FluidBlock& b) const
  {
    const bool touchW = fadeLen[0] >= b.min_pos[0];
    const bool touchE = fadeLen[0] >= ext[0] - b.max_pos[0];
    const bool touchS = fadeLen[1] >= b.min_pos[1];
    const bool touchN = fadeLen[1] >= ext[1] - b.max_pos[1];
    const bool touchB = fadeLen[2] >= b.min_pos[2];
    const bool touchF = fadeLen[2] >= ext[2] - b.max_pos[2];
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
  KernelFadeOut(const Real buf[3], const Real extent[3]) :
  ext{extent[0],extent[1],extent[2]}, fadeLen{buf[0],buf[1],buf[2]},
  iFade{1/(buf[0]+EPS), 1/(buf[1]+EPS), 1/(buf[2]+EPS)} {}

  void operator()(const BlockInfo& info, FluidBlock& b) const
  {
    if( _is_touching(b) )
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const Real FADE = fade(info, ix, iy, iz);
      b(ix,iy,iz).u *= FADE; b(ix,iy,iz).v *= FADE; b(ix,iy,iz).w *= FADE;
    }
  }
};

}

void FadeOut::operator()(const double dt)
{
  sim.startProfiler("FadeOut Kernel");
  #pragma omp parallel
  {
    KernelFadeOut kernel(sim.fadeOutLengthU, sim.extent.data());
    #pragma omp for schedule(static)
    for (size_t i=0; i<vInfo.size(); i++)
      kernel(vInfo[i], *(FluidBlock*) vInfo[i].ptrBlock);
  }
  sim.stopProfiler();
  check("FadeOut");
}

CubismUP_3D_NAMESPACE_END
