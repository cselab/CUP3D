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

#if 0
namespace EllipsoidObstacle
{
// code from http://www.geometrictools.com/
//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1 with e0 >= e1 >= e2.
// The query point is (y0,y1,y2) with y0 >= 0, y1 >= 0, and y2 >= 0.  The
// function returns the distance from the query point to the ellipsoid.  It
// also computes the ellipsoid point (x0,x1,x2) in the first octant that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------

static Real DistancePointEllipsoidSpecial(const Real e[3],const Real y[3],Real x[3])
{
  if (y[2] > (Real)0)
  {
    if (y[1] > (Real)0)
    {
      if (y[0] > (Real)0)
      {
        // Bisect to compute the root of F(t) for t >= -e2*e2.
        const Real esqr[3] = { e[0]*e[0], e[1]*e[1], e[2]*e[2] };
        const Real ey[3] = { e[0]*y[0], e[1]*y[1], e[2]*y[2] };
        Real t0 = -esqr[2] + ey[2];
        Real t1 = -esqr[2] + std::sqrt(ey[0]*ey[0] +ey[1]*ey[1] +ey[2]*ey[2]);
        Real t = t0;
        const int imax = 2*std::numeric_limits<Real>::max_exponent;
        for (int i = 0; i < imax; ++i) {
          t = (t0 + t1)/2;
          static constexpr Real eps = std::numeric_limits<Real>::epsilon();
          if (std::fabs(t-t0)<eps || std::fabs(t-t1)<eps) break;
          const Real r[3] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]),
              ey[2]/(t + esqr[2]) };
          const Real f = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1;
          if (f > 0) t0 = t;
          else if (f < 0) t1 = t;
          else break;
        }

        x[0] = esqr[0]*y[0]/(t + esqr[0]);
        x[1] = esqr[1]*y[1]/(t + esqr[1]);
        x[2] = esqr[2]*y[2]/(t + esqr[2]);
        const Real d[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
        return std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      }
      else  // y0 == 0
      {
        x[0] = 0;
        const Real etmp[2] = { e[1], e[2] };
        const Real ytmp[2] = { y[1], y[2] };
        Real xtmp[2];
        const Real distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
        x[1] = xtmp[0];
        x[2] = xtmp[1];
        return distance;
      }
    }
    else  // y1 == 0
    {
      x[1] = (Real)0;
      if (y[0] > (Real)0)
      {
        const Real etmp[2] = { e[0], e[2] };
        const Real ytmp[2] = { y[0], y[2] };
        const Real xtmp[2];
        const Real distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
        x[0] = xtmp[0];
        x[2] = xtmp[1];
        return distance;
      }
      else  // y0 == 0
      {
        x[0] = (Real)0;
        x[2] = e[2];
        return std::fabs(y[2] - e[2]);
      }
    }
  }
  else  // y2 == 0
  {
    const Real denom[2] = { e[0]*e[0] - e[2]*e[2], e[1]*e[1] - e[2]*e[2] };
    const Real ey[2] = { e[0]*y[0], e[1]*y[1] };
    if (ey[0] < denom[0] && ey[1] < denom[1])
    {
      // (y0,y1) is inside the axis-aligned bounding rectangle of the
      // subellipse.  This intermediate test is designed to guard
      // against the division by zero when e0 == e2 or e1 == e2.
      const Real xde[2] = { ey[0]/denom[0], ey[1]/denom[1] };
      const Real xdesqr[2] = { xde[0]*xde[0], xde[1]*xde[1] };
      const Real discr = 1 - xdesqr[0] - xdesqr[1];
      if (discr > (Real)0)
      {
        // (y0,y1) is inside the subellipse.  The closest ellipsoid
        // point has x2 > 0.
        x[0] = e[0]*xde[0];
        x[1] = e[1]*xde[1];
        x[2] = e[2]*sqrt(discr);
        const Real d[2] = { x[0] - y[0], x[1] - y[1] };
        return std::sqrt(d[0]*d[0] + d[1]*d[1] + x[2]*x[2]);
      } else {
        // (y0,y1) is outside the subellipse.  The closest ellipsoid
        // point has x2 == 0 and is on the domain-boundary ellipse
        // (x0/e0)^2 + (x1/e1)^2 = 1.
        x[2] = (Real)0;
        return DistancePointEllipseSpecial(e, y, x);
      }
    } else {
      // (y0,y1) is outside the subellipse.  The closest ellipsoid
      // point has x2 == 0 and is on the domain-boundary ellipse
      // (x0/e0)^2 + (x1/e1)^2 = 1.
      x[2] = (Real)0;
      return DistancePointEllipseSpecial(e, y, x);
    }
  }
}
//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1.  The query point is
// (y0,y1,y2).  The function returns the distance from the query point to the
// ellipsoid.   It also computes the ellipsoid point (x0,x1,x2) that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------
static Real DistancePointEllipsoid (const Real e[3], const Real y[3], Real x[3])
{
  // Determine reflections for y to the first octant.
  const bool reflect[3] = {y[0]<0, y[1]<0, y[2]<0};

  // Determine the axis order for decreasing extents.
  int permute[3];
  if (e[0] < e[1]) {
    if (e[2] < e[0]) {
      permute[0] = 1;  permute[1] = 0;  permute[2] = 2;
    }
    else if (e[2] < e[1]) {
      permute[0] = 1;  permute[1] = 2;  permute[2] = 0;
    }
    else {
      permute[0] = 2;  permute[1] = 1;  permute[2] = 0;
    }
  } else {
    if (e[2] < e[1]) {
      permute[0] = 0;  permute[1] = 1;  permute[2] = 2;
    }
    else if (e[2] < e[0]) {
      permute[0] = 0;  permute[1] = 2;  permute[2] = 1;
    }
    else {
      permute[0] = 2;  permute[1] = 0;  permute[2] = 1;
    }
  }

  int invpermute[3];
  for (int i = 0; i < 3; ++i) invpermute[permute[i]] = i;

  Real locE[3], locY[3];
  for (i = 0; i < 3; ++i) {
    j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if (reflect[j]) locY[i] = -locY[i];
  }

  Real locX[3];
  const Real distance = DistancePointEllipsoidSpecial(locE, locY, locX);

  // Restore the axis order and reflections.
  for (i = 0; i < 3; ++i) {
    j = invpermute[i];
    if (reflect[j]) locX[j] = -locX[j];
    x[i] = locX[j];
  }

  return distance;
}

struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real e0,e1,e2,h;
  const Real maxAxis = std::max(std::max(e0,e1),e2);
  const Real position[3], quaternion[4];
  Real box[3][2] = {
    {position[0] - 2*(maxAxis + 2*h), position[0] + 2*(maxAxis + 2*h)},
    {position[1] - 2*(maxAxis + 2*h), position[1] + 2*(maxAxis + 2*h)},
    {position[2] - 2*(maxAxis + 2*h), position[2] + 2*(maxAxis + 2*h)}
  };
  const Real w = quaternion[0], x = quaternion[1], y = quaternion[2], z = quaternion[3];
  const Real Rmatrix[3][3] = {
      {1-2*(y*y+z*z),  2*(x*y+z*w),    2*(x*z-y*w)},
      {2*(x*y-z*w),    1-2*(x*x+z*z),  2*(y*z+x*w)},
      {2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x*x+y*y)}
  };

  FillBlocks(const Real _e0, const Real _e1, const Real _e2, const Real _h,
    const Real p[3], const Real q[4]) : e0(_e0), e1(_e1), e2(_e2), h(_h),
    position{p[0],p[1],p[2]}, quaternion{q[0],q[1],q[2],q[3]} { }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    Real intersection[3][2] = {
      std::max(min_pos[0], box[0][0]), std::min(max_pos[0], box[0][1]),
      std::max(min_pos[1], box[1][0]), std::min(max_pos[1], box[1][1]),
      std::max(min_pos[2], box[2][0]), std::min(max_pos[2], box[2][1])
    };
    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool isTouching(const BlockInfo& info, const int buffer_dx=0) const
  {
    Real min_pos[3], max_pos[3];
    info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i) {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return _is_touching(min_pos,max_pos);
  }

  inline Real signedDistance(const Real xo, const Real yo, const Real zo) const
  {
    const Real p[3] = {xo-position[0], yo-position[1], zo-position[2]};
    // rotate
    const Real t[3] = {
        Rmatrix[0][0]*p[0] + Rmatrix[0][1]*p[1] + Rmatrix[0][2]*p[2],
        Rmatrix[1][0]*p[0] + Rmatrix[1][1]*p[1] + Rmatrix[1][2]*p[2],
        Rmatrix[2][0]*p[0] + Rmatrix[2][1]*p[1] + Rmatrix[2][2]*p[2]
    };
    // find distance
    const Real e[3] = {e0, e1, e2};
    Real xs[3];
    const Real dist = DistancePointEllipsoid (e, t, xs);
    const Real Dcentre = t[0]*t[0]+t[1]*t[1]+t[2]*t[2];
    const Real Dsurf = xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2];
    const int sign = Dcentre > Dsurf ? 1 : -1;
  }
};
}
#endif
#endif
