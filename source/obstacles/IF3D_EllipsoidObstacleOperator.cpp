//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "obstacles/IF3D_EllipsoidObstacleOperator.h"
#include "obstacles/extra/IF3D_ObstacleLibrary.h"
#include "Cubism/ArgumentParser.h"

namespace EllipsoidObstacle
{
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
static constexpr int IMAX = 2*std::numeric_limits<Real>::max_exponent;

static Real distPointEllipseSpecial(const Real e[2], const Real y[2], Real x[2])
{
  if (y[1] > 0) {
    if (y[0] > 0) {
      // Bisect to compute the root of F(t) for t >= -e1*e1.
      const Real esqr[2] = { e[0]*e[0], e[1]*e[1] };
      const Real ey[2] = { e[0]*y[0], e[1]*y[1] };
      Real t0 = -esqr[1] + ey[1];
      Real t1 = -esqr[1] + std::sqrt(ey[0]*ey[0] + ey[1]*ey[1]);
      Real t = t0;
      for (int i = 0; i < IMAX; ++i) {
        t = 0.5 * (t0 + t1);
        if (std::fabs(t-t0)<EPS || std::fabs(t-t1)<EPS) break;
        const Real r[2] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]) };
        const Real f = r[0]*r[0] + r[1]*r[1] - 1;
        if (f > 0) t0 = t;
        else if (f < 0) t1 = t;
        else break;
      }
      x[0] = esqr[0]*y[0]/(t + esqr[0]);
      x[1] = esqr[1]*y[1]/(t + esqr[1]);
      const Real d[2] = { x[0] - y[0], x[1] - y[1] };
      return std::sqrt(d[0]*d[0] + d[1]*d[1]);
    } else { // y0 == 0
      x[0] = 0;
      x[1] = e[1];
      return std::fabs(y[1] - e[1]);
    }
  } else { // y1 == 0
    const Real denom0 = e[0]*e[0] - e[1]*e[1];
    const Real e0y0 = e[0]*y[0];
    if (e0y0 < denom0) {
      // y0 is inside the subinterval.
      const Real x0de0 = e0y0/denom0;
      const Real x0de0sqr = x0de0*x0de0;
      x[0] = e[0]*x0de0;
      x[1] = e[1]*sqrt(fabs((Real)1 - x0de0sqr));
      const Real d0 = x[0] - y[0];
      return std::sqrt(d0*d0 + x[1]*x[1]);
    } else {
      // y0 is outside the subinterval.  The closest ellipse point has
      // x1 == 0 and is on the domain-boundary interval (x0/e0)^2 = 1.
      x[0] = e[0];
      x[1] = 0;
      return std::fabs(y[0] - e[0]);
    }
  }
}

//----------------------------------------------------------------------------
// The ellipse is (x0/e0)^2 + (x1/e1)^2 = 1.  The query point is (y0,y1).
// The function returns the distance from the query point to the ellipse.
// It also computes the ellipse point (x0,x1) that is closest to (y0,y1).
//----------------------------------------------------------------------------
/*
static Real distancePointEllipse (const Real e[2], const Real y[2], Real x[2])
{
  // Determine reflections for y to the first quadrant.
  const bool reflect[2] = {y[0] < 0, y[1] < 0};
  // Determine the axis order for decreasing extents.
  const int permute[2] = {e[0]<e[1] ? 1 : 0, e[0]<e[1] ? 0 : 1};
  int invpermute[2];
  for (int i = 0; i < 2; ++i) invpermute[permute[i]] = i;
  Real locE[2], locY[2];
  for (int i = 0; i < 2; ++i) {
    const int j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if (reflect[j]) locY[i] = -locY[i];
  }

  Real locX[2];
  const Real distance = distPointEllipseSpecial(locE, locY, locX);

  // Restore the axis order and reflections.
  for (int i = 0; i < 2; ++i) {
    const int j = invpermute[i];
    if (reflect[j]) locX[j] = -locX[j];
    x[i] = locX[j];
  }
  return distance;
}
*/
// code from http://www.geometrictools.com/
//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1 with e0 >= e1 >= e2.
// The query point is (y0,y1,y2) with y0 >= 0, y1 >= 0, and y2 >= 0.  The
// function returns the distance from the query point to the ellipsoid.  It
// also computes the ellipsoid point (x0,x1,x2) in the first octant that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------

static Real distPointEllipsoidSpecial(const Real e[3],const Real y[3],Real x[3])
{
  if (y[2] > 0) {
    if (y[1] > 0) {
      if (y[0] > 0) {
        // Bisect to compute the root of F(t) for t >= -e2*e2.
        const Real esq[3] = { e[0]*e[0], e[1]*e[1], e[2]*e[2] };
        const Real ey[3] =  { e[0]*y[0], e[1]*y[1], e[2]*y[2] };
        Real t0 = -esq[2] + ey[2];
        Real t1 = -esq[2] + std::sqrt(ey[0]*ey[0] +ey[1]*ey[1] +ey[2]*ey[2]);
        Real t = t0;
        for (int i = 0; i < IMAX; ++i) {
          t = 0.5 * (t0 + t1);
          if (std::fabs(t-t0)<EPS || std::fabs(t-t1)<EPS) break;
          const Real r[3]= {ey[0]/(t+esq[0]),ey[1]/(t+esq[1]),ey[2]/(t+esq[2])};
          const Real f = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - 1;
          if (f > 0) t0 = t;
          else if (f < 0) t1 = t;
          else break;
        }
        x[0] = esq[0]*y[0]/(t + esq[0]);
        x[1] = esq[1]*y[1]/(t + esq[1]);
        x[2] = esq[2]*y[2]/(t + esq[2]);
        const Real d[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
        return std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      } else { // y0 == 0
        x[0] = 0;
        const Real etmp[2] = { e[1], e[2] };
        const Real ytmp[2] = { y[1], y[2] };
        Real xtmp[2];
        const Real distance = distPointEllipseSpecial(etmp, ytmp, xtmp);
        x[1] = xtmp[0];
        x[2] = xtmp[1];
        return distance;
      }
    } else { // y1 == 0
      x[1] = 0;
      if (y[0] > 0) {
        const Real etmp[2] = { e[0], e[2] };
        const Real ytmp[2] = { y[0], y[2] };
        Real xtmp[2];
        const Real distance = distPointEllipseSpecial(etmp, ytmp, xtmp);
        x[0] = xtmp[0];
        x[2] = xtmp[1];
        return distance;
      } else { // y0 == 0
        x[0] = 0;
        x[2] = e[2];
        return std::fabs(y[2] - e[2]);
      }
    }
  } else { // y2 == 0
    const Real denom[2] = { e[0]*e[0] - e[2]*e[2], e[1]*e[1] - e[2]*e[2] };
    const Real ey[2] = { e[0]*y[0], e[1]*y[1] };
    if (ey[0] < denom[0] && ey[1] < denom[1]) {
      // (y0,y1) is inside the axis-aligned bounding rectangle of the
      // subellipse.  This intermediate test is designed to guard
      // against the division by zero when e0 == e2 or e1 == e2.
      const Real xde[2] = { ey[0]/denom[0], ey[1]/denom[1] };
      const Real xdesqr[2] = { xde[0]*xde[0], xde[1]*xde[1] };
      const Real discr = 1 - xdesqr[0] - xdesqr[1];
      if (discr > 0) {
        // (y0,y1) is inside the subellipse.  The closest ellipsoid
        // point has x2 > 0.
        x[0] = e[0]*xde[0];
        x[1] = e[1]*xde[1];
        x[2] = e[2]*std::sqrt(discr);
        const Real d[2] = { x[0] - y[0], x[1] - y[1] };
        return std::sqrt(d[0]*d[0] + d[1]*d[1] + x[2]*x[2]);
      } else {
        // (y0,y1) is outside the subellipse.  The closest ellipsoid
        // point has x2 == 0 and is on the domain-boundary ellipse
        // (x0/e0)^2 + (x1/e1)^2 = 1.
        x[2] = 0;
        return distPointEllipseSpecial(e, y, x);
      }
    } else {
      // (y0,y1) is outside the subellipse.  The closest ellipsoid
      // point has x2 == 0 and is on the domain-boundary ellipse
      // (x0/e0)^2 + (x1/e1)^2 = 1.
      x[2] = 0;
      return distPointEllipseSpecial(e, y, x);
    }
  }
}

//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1.  The query point is
// (y0,y1,y2).  The function returns the distance from the query point to the
// ellipsoid.   It also computes the ellipsoid point (x0,x1,x2) that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------
static Real distancePointEllipsoid(const Real e[3], const Real y[3], Real x[3])
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
  for (int i = 0; i < 3; ++i) {
    const int j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if (reflect[j]) locY[i] = -locY[i];
  }

  Real locX[3];
  const Real distance = distPointEllipsoidSpecial(locE, locY, locX);

  // Restore the axis order and reflections.
  for (int i = 0; i < 3; ++i) {
    const int j = invpermute[i];
    if (reflect[j]) locX[j] = -locX[j];
    x[i] = locX[j];
  }

  return distance;
}

struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real e0, e1, e2, h;
  const Real maxAxis = std::max({e0, e1, e2});
  const Real position[3], quaternion[4];
  const Real box[3][2] = {
    {position[0] - 2*(maxAxis + 2*h), position[0] + 2*(maxAxis + 2*h)},
    {position[1] - 2*(maxAxis + 2*h), position[1] + 2*(maxAxis + 2*h)},
    {position[2] - 2*(maxAxis + 2*h), position[2] + 2*(maxAxis + 2*h)}
  };
  const Real w=quaternion[0], x=quaternion[1], y=quaternion[2], z=quaternion[3];
  const Real Rmatrix[3][3] = {
      {1-2*(y*y+z*z),   2*(x*y+z*w),   2*(x*z-y*w)},
      {  2*(x*y-z*w), 1-2*(x*x+z*z),   2*(y*z+x*w)},
      {  2*(x*z+y*w),   2*(y*z-x*w), 1-2*(x*x+y*y)}
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
    const Real dist = distancePointEllipsoid (e, t, xs);
    const Real Dcentre = t[0]*t[0]+t[1]*t[1]+t[2]*t[2];
    const Real Dsurf = xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2];
    const int sign = Dcentre > Dsurf ? 1 : -1;
    return dist * sign;
  }
};
}

IF3D_EllipsoidObstacleOperator::IF3D_EllipsoidObstacleOperator(
    SimulationData& s, ArgumentParser& p )
    : IF3D_ObstacleOperator(s, p), radius(0.5*length)
{
  e0 = p("-aspectRatioX").asDouble(1) * radius;
  e1 = p("-aspectRatioY").asDouble(1) * radius;
  e2 = p("-aspectRatioZ").asDouble(1) * radius;
  accel_decel = p("-accel").asBool(false);
  if(accel_decel) {
    if(not bForcedInSimFrame[0]) {
      printf("Warning: sphere was not set to be forced in x-dir, yet the accel_decel pattern is active.\n");
    }
    umax = p("-xvel").asDouble(0.0);
    tmax = p("-T").asDouble(1.);
  }
}

void IF3D_EllipsoidObstacleOperator::create()
{
  const Real h = vInfo[0].h_gridpoint;
  const EllipsoidObstacle::FillBlocks K(e0,e1,e2, h, position, quaternion);

  create_base<EllipsoidObstacle::FillBlocks>(K);
}

void IF3D_EllipsoidObstacleOperator::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}


void IF3D_EllipsoidObstacleOperator::computeVelocities()
{
  IF3D_ObstacleOperator::computeVelocities();

  if(accel_decel) {
    if(sim.time<tmax) transVel[0] = umax*sim.time/tmax;
    else if (sim.time<2*tmax) transVel[0] = umax*(2*tmax-sim.time)/tmax;
    else transVel[0] = 0;
  }
}
