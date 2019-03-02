//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef IncompressibleFluids3D_IF3D_ObstacleLibrary_h
#define IncompressibleFluids3D_IF3D_ObstacleLibrary_h

#include "Definitions.h"
#include "obstacles/extra/IF2D_Interpolation1D.h"

/*
 * TODO: Refactor `SphereObstacle` and `DCylinderObstacle` to derive from
 *       `FillBlocksBase`. Also, move them to corresponding obstacle .cpp files.
 */

/*
 * A base class for FillBlocks classes.
 *
 * Derived classes should be implemented as (*):
 *      class FillBlocksFOO : FillBlocksBase<FillBlocksFOO> {
 *          (...)
 *      };
 *
 * and are required to implement following member functions:
 *
 * bool isTouching(const BlockInfo &info, int buffer_dx = 0) const;
 *      Returns if the given blocks intersects or touches the object.
 *      False positives are acceptable.
 *
 * Real signedDistance(Real x, Real y, Real z) const;
 *      Returns the signed distance of the given point from the surface of the
 *      object. Positive number stands for inside, negative for outside.
 *
 *
 * (*) That way the base class is able to access functions of the derived class
 *     in compile-time (without using virtual functions). For more info, see:
 *     https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 */
template <typename Derived>
struct FillBlocksBase
{
using CHIMAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX];
#define DERIVED (static_cast<const Derived *>(this))
    void operator()(const BlockInfo &info, ObstacleBlock* const oBlock) const
    {
      // TODO: Remove `isTouching` check and verify that all dependencies are
      //       using this function properly.
      if (!DERIVED->isTouching(info)) return;
      CHIMAT & __restrict__ SDF = oBlock->sdf;
      FluidBlock &b = *(FluidBlock *)info.ptrBlock;
      for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
      for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
      for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
        Real p[3];
        info.pos(p, ix, iy, iz);
        const Real dist = DERIVED->signedDistance(p[0], p[1], p[2]);
        SDF[iz][iy][ix] = dist;
        // negative outside of the obstacle, therefore max = minimal distance.
        b(ix,iy,iz).tmpU = std::max(dist, b(ix,iy,iz).tmpU);
      }
    }
#undef DERIVED
};

namespace TorusObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real big_r, small_r, h, position[3];
  Real sphere_box[3][2] = {
    {position[0] - 2*(small_r + 2*h),    position[0] + 2*(small_r + 2*h)},
    {position[1] -2*(big_r+small_r+2*h), position[1] +2*(big_r+small_r+2*h)},
    {position[2] -2*(big_r+small_r+2*h), position[2] +2*(big_r+small_r+2*h)}
  };

  FillBlocks(Real _big_r, Real _small_r, Real _h,  Real p[3]) :
    big_r(_big_r), small_r(_small_r), h(_h), position{p[0],p[1],p[2]} { }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    using std::max; using std::min;
    const Real intersection[3][2] = {
        max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
        max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
        max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool isTouching(const BlockInfo& info, const int buffer_dx=0) const
  {
    Real min_pos[3], max_pos[3]; info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i) {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return _is_touching(min_pos,max_pos);
  }

  inline Real signedDistance(const Real xo, const Real yo, const Real zo) const
  {
    const Real t[3] = { xo - position[0], yo - position[1], zo - position[2] };
    const Real r = std::sqrt(t[1]*t[1] + t[2]*t[2]);
    if (r > 0) {
      const Real c[3] = { 0, big_r*t[1]/r, big_r*t[2]/r };
      const Real d = std::pow(t[0]-c[0],2) + std::pow(t[1]-c[1],2) + std::pow(t[2]-c[2],2);
      return std::sqrt(d) - small_r;
    }
    else return -1; // very far??? no else in original code: chi = 0
  }
};
}
#endif
