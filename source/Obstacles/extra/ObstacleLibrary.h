//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef CubismUP_3D_ObstacleLibrary_h
#define CubismUP_3D_ObstacleLibrary_h

#include "../../Definitions.h"

CubismUP_3D_NAMESPACE_BEGIN

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
 * bool isTouching(const BlockInfo &info, const FluidBlock &block) const;
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
  void operator()(const cubism::BlockInfo &info, ObstacleBlock* const o) const
  {
    // TODO: Remove `isTouching` check and verify that all dependencies are
    //       using this function properly.
    FluidBlock &b = *(FluidBlock *)info.ptrBlock;
    if (!derived()->isTouching(info, b)) return;
    auto & SDFLAB = o->sdfLab;
    for (int iz = -1; iz < FluidBlock::sizeZ+1; ++iz)
    for (int iy = -1; iy < FluidBlock::sizeY+1; ++iy)
    for (int ix = -1; ix < FluidBlock::sizeX+1; ++ix) {
      Real p[3];
      info.pos(p, ix, iy, iz);
      const Real dist = derived()->signedDistance(p[0], p[1], p[2]);
      SDFLAB[iz+1][iy+1][ix+1] = dist;
    }
  }

private:
  const Derived* derived() const noexcept
  {
    return static_cast<const Derived *>(this);
  }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ObstacleLibrary_h
