//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef IncompressibleFluids3D_IF3D_ObstacleLibrary_h
#define IncompressibleFluids3D_IF3D_ObstacleLibrary_h
#include "IF2D_Interpolation1D.h"
#include "GenericOperator.h"


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
struct FillBlocksBase {
#define DERIVED (static_cast<const Derived *>(this))
    void operator()(const BlockInfo &info, ObstacleBlock * const o) const
    {
        /*
         * Based on:
         *
         * Towers, John D.
         * "Finite difference methods for approximating Heaviside functions."
         * Journal of Computational Physics 228.9 (2009): 3478-3489.
         *
         * https://www.sciencedirect.com/science/article/pii/S0021999109000576
         */

        // TODO: Remove `isTouching` check and verify that all dependencies are
        //       using this function properly.
        if (!DERIVED->isTouching(info))
            return;

        FluidBlock &b = *(FluidBlock *)info.ptrBlock;
        const Real h = info.h_gridpoint;
        const Real inv2h = (Real)0.5 / h;
        const Real fac1 = (Real)0.5 * h * h;
        static constexpr Real eps = std::numeric_limits<Real>::epsilon();

        for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
        for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
        for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
            Real p[3];
            info.pos(p, ix, iy, iz);
            const Real dist = DERIVED->signedDistance(p[0], p[1], p[2]);

            if (dist > 2 * h || dist < -2 * h) { //2 should be safe
                const Real H = dist > 0 ? 1 : 0;
                b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
                o->write(ix, iy, iz, H, 0, 0, 0, 0, 0);
                continue;
            }

            // This should be cached in a matrix, recalc might be expensive.
            const Real distPx = DERIVED->signedDistance(p[0] + h, p[1], p[2]);
            const Real distMx = DERIVED->signedDistance(p[0] - h, p[1], p[2]);
            const Real distPy = DERIVED->signedDistance(p[0], p[1] + h, p[2]);
            const Real distMy = DERIVED->signedDistance(p[0], p[1] - h, p[2]);
            const Real distPz = DERIVED->signedDistance(p[0], p[1], p[2] + h);
            const Real distMz = DERIVED->signedDistance(p[0], p[1], p[2] - h);

            const Real IplusX = distPx < 0 ? 0 : distPx;
            const Real IminuX = distMx < 0 ? 0 : distMx;
            const Real IplusY = distPy < 0 ? 0 : distPy;
            const Real IminuY = distMy < 0 ? 0 : distMy;
            const Real IplusZ = distPz < 0 ? 0 : distPz;
            const Real IminuZ = distMz < 0 ? 0 : distMz;

            const Real HplusX = std::fabs(distPx)<eps ? (Real)0.5 :
                                                        (distPx < 0 ? 0 : 1);
            const Real HminuX = std::fabs(distMx)<eps ? (Real)0.5 :
                                                        (distMx < 0 ? 0 : 1);
            const Real HplusY = std::fabs(distPy)<eps ? (Real)0.5 :
                                                        (distPy < 0 ? 0 : 1);
            const Real HminuY = std::fabs(distMy)<eps ? (Real)0.5 :
                                                        (distMy < 0 ? 0 : 1);
            const Real HplusZ = std::fabs(distPz)<eps ? (Real)0.5 :
                                                        (distPz < 0 ? 0 : 1);
            const Real HminuZ = std::fabs(distMz)<eps ? (Real)0.5 :
                                                        (distMz < 0 ? 0 : 1);

            // all would be multiplied by 0.5/h, simplifies out later
            const Real gradUX = inv2h * (distPx - distMx);
            const Real gradUY = inv2h * (distPy - distMy);
            const Real gradUZ = inv2h * (distPz - distMz);
            const Real gradIX = inv2h * (IplusX - IminuX);
            const Real gradIY = inv2h * (IplusY - IminuY);
            const Real gradIZ = inv2h * (IplusZ - IminuZ);
            const Real gradHX = HplusX - HminuX;
            const Real gradHY = HplusY - HminuY;
            const Real gradHZ = HplusZ - HminuZ;

            const Real gradUSq = gradUX * gradUX + gradUY * gradUY + gradUZ * gradUZ + eps;
            const Real numH    = gradIX * gradUX + gradIY * gradUY + gradIZ * gradUZ;
            const Real numD    = gradHX * gradUX + gradHY * gradUY + gradHZ * gradUZ;
            const Real Delta   = numD / gradUSq;
            const Real H       = numH / gradUSq;

            o->write(ix, iy, iz, H, Delta, gradUX, gradUY, gradUZ, fac1);
            b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
#ifndef NDEBUG
            if (H < 0 || H > 1 + eps)
                printf("invalid H?: %f %f %f: %10.10e\n", p[0], p[1], p[2], H);
#endif
        }
    }
#undef DERIVED
};

namespace SphereObstacle
{
struct FillBlocks
{
  const Real radius,safe_radius;
  const double sphere_position[3];
  Real sphere_box[3][2];

  void _find_sphere_box()
  {
    sphere_box[0][0] = sphere_position[0] - safe_radius;
    sphere_box[0][1] = sphere_position[0] + safe_radius;
    sphere_box[1][0] = sphere_position[1] - safe_radius;
    sphere_box[1][1] = sphere_position[1] + safe_radius;
    sphere_box[2][0] = sphere_position[2] - safe_radius;
    sphere_box[2][1] = sphere_position[2] + safe_radius;
  }

  FillBlocks(const Real _radius, const Real max_dx, const double pos[3]):
    radius(_radius), safe_radius(radius+4*max_dx), sphere_position{pos[0],pos[1],pos[2]}
  {
    _find_sphere_box();
  }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    using std::min;
    using std::max;
    Real intersection[3][2] = {
        max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
        max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
        max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool _is_touching(const BlockInfo& info, const int buffer_dx = 0) const
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

  inline Real distanceToSphere(const Real x, const Real y, const Real z) const
  {
    return radius - std::sqrt(x*x+y*y+z*z); // pos inside, neg outside
  }

  inline Real sign(const Real& val) const {
    return (0. < val) - (val < 0.);
  }

  void operator()(const BlockInfo&info, FluidBlock&b, ObstacleBlock*const o) const
  {
    if(_is_touching(info)) {
      const Real h = info.h_gridpoint;
      const Real inv2h = 0.5/h;
      const Real fac1 = 0.5*h*h;
      static constexpr Real eps = std::numeric_limits<Real>::epsilon();

      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
      for(int iy=0; iy<FluidBlock::sizeY; iy++)
      for(int ix=0; ix<FluidBlock::sizeX; ix++) {
        Real p[3];
        info.pos(p, ix, iy, iz);
        const Real x = p[0]-sphere_position[0];
        const Real y = p[1]-sphere_position[1];
        const Real z = p[2]-sphere_position[2];
        const Real dist = distanceToSphere(x,y,z);

        if(dist > 2*h || dist < -2*h) { //2 should be safe
          const Real H = dist > 0 ? 1.0 : 0.0;
          b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
          o->write(ix,iy,iz,H,0,0,0,0,0);
          continue;
        }

        const Real distPx = distanceToSphere(x+h,y,z);
        const Real distMx = distanceToSphere(x-h,y,z);
        const Real distPy = distanceToSphere(x,y+h,z);
        const Real distMy = distanceToSphere(x,y-h,z);
        const Real distPz = distanceToSphere(x,y,z+h);
        const Real distMz = distanceToSphere(x,y,z-h);
        const Real IplusX = distPx < 0 ? 0 : distPx;
        const Real IminuX = distMx < 0 ? 0 : distMx;
        const Real IplusY = distPy < 0 ? 0 : distPy;
        const Real IminuY = distMy < 0 ? 0 : distMy;
        const Real IplusZ = distPz < 0 ? 0 : distPz;
        const Real IminuZ = distMz < 0 ? 0 : distMz;

        const Real HplusX = std::fabs(distPx)<eps ? (Real)0.5 :
                                                    (distPx < 0 ? 0 : 1);
        const Real HminuX = std::fabs(distMx)<eps ? (Real)0.5 :
                                                    (distMx < 0 ? 0 : 1);
        const Real HplusY = std::fabs(distPy)<eps ? (Real)0.5 :
                                                    (distPy < 0 ? 0 : 1);
        const Real HminuY = std::fabs(distMy)<eps ? (Real)0.5 :
                                                    (distMy < 0 ? 0 : 1);
        const Real HplusZ = std::fabs(distPz)<eps ? (Real)0.5 :
                                                    (distPz < 0 ? 0 : 1);
        const Real HminuZ = std::fabs(distMz)<eps ? (Real)0.5 :
                                                    (distMz < 0 ? 0 : 1);

        // all would be multiplied by 0.5/h, simplifies out later
        const Real gradUX = inv2h*(distPx - distMx);
        const Real gradUY = inv2h*(distPy - distMy);
        const Real gradUZ = inv2h*(distPz - distMz);
        const Real gradIX = inv2h*(IplusX - IminuX);
        const Real gradIY = inv2h*(IplusY - IminuY);
        const Real gradIZ = inv2h*(IplusZ - IminuZ);
        const Real gradHX = (HplusX - HminuX);
        const Real gradHY = (HplusY - HminuY);
        const Real gradHZ = (HplusZ - HminuZ);

        const Real gradUSq = gradUX*gradUX +gradUY*gradUY +gradUZ*gradUZ +eps;
        const Real numH    = gradIX*gradUX +gradIY*gradUY +gradIZ*gradUZ;
        const Real numD    = gradHX*gradUX +gradHY*gradUY +gradHZ*gradUZ;
        const Real Delta   = numD / gradUSq;
        const Real H       = numH / gradUSq;

        o->write(ix, iy, iz, H, Delta, gradUX, gradUY, gradUZ, fac1);
        b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
        #ifndef NDEBUG
        if(H<0||H>1+eps) printf("invalid H?: %f %f %f: %10.10e\n",x,y,z,H);
        #endif
      }
    }
  }
};
}

namespace DCylinderObstacle
{
struct FillBlocks
{
  const Real radius, halflength, safety;
  const double position[3];
  Real box[3][2];

  FillBlocks(const Real r, const Real halfl, const Real h, const double p[3]):
  radius(r), halflength(halfl), safety(2*h), position{p[0],p[1],p[2]}
  {
    box[0][0] = position[0] - radius     - safety;
    box[0][1] = position[0]              + safety;
    box[1][0] = position[1] - radius     - safety;
    box[1][1] = position[1] + radius     + safety;
    box[2][0] = position[2] - halflength - safety;
    box[2][1] = position[2] + halflength + safety;
  }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    using std::min;
    using std::max;
    Real intersection[3][2] = {
        max(min_pos[0], box[0][0]), min(max_pos[0], box[0][1]),
        max(min_pos[1], box[1][0]), min(max_pos[1], box[1][1]),
        max(min_pos[2], box[2][0]), min(max_pos[2], box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool _is_touching(const BlockInfo& info, const int buffer_dx = 0) const
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

  inline Real distance(const Real x, const Real y, const Real z) const
  {
    const Real planeDist = std::min( -x, radius-std::sqrt(x*x+y*y) );
    const Real vertiDist = halflength - std::fabs(z);
    return std::min(planeDist, vertiDist);
  }

  void operator()(const BlockInfo&info, FluidBlock&b, ObstacleBlock*const o) const
  {
    if(_is_touching(info)) {
      const double h=info.h_gridpoint;
      const Real eps = std::numeric_limits<Real>::epsilon();
      const Real inv2h = 0.5/h, fac1 = 0.5*h*h;

      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
      for(int iy=0; iy<FluidBlock::sizeY; iy++)
      for(int ix=0; ix<FluidBlock::sizeX; ix++) {
        Real p[3];
        info.pos(p, ix, iy, iz);
        const Real x=p[0]-position[0], y=p[1]-position[1], z=p[2]-position[2];
        const Real dist = distance(x,y,z);
        if(dist > 2*h || dist < -2*h) { //2 should be safe
          const Real H = dist > 0 ? 1.0 : 0.0;
          b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
          o->write(ix,iy,iz,H,0,0,0,0,0);
          continue;
        }

        const Real distPx = distance(x+h,y,z);
        const Real distMx = distance(x-h,y,z);
        const Real distPy = distance(x,y+h,z);
        const Real distMy = distance(x,y-h,z);
        const Real distPz = distance(x,y,z+h);
        const Real distMz = distance(x,y,z-h);
        const Real IplusX = distPx < 0 ? 0 : distPx;
        const Real IminuX = distMx < 0 ? 0 : distMx;
        const Real IplusY = distPy < 0 ? 0 : distPy;
        const Real IminuY = distMy < 0 ? 0 : distMy;
        const Real IplusZ = distPz < 0 ? 0 : distPz;
        const Real IminuZ = distMz < 0 ? 0 : distMz;

        const Real HplusX = std::fabs(distPx)<eps ? (Real)0.5 :
                                                    (distPx < 0 ? 0 : 1);
        const Real HminuX = std::fabs(distMx)<eps ? (Real)0.5 :
                                                    (distMx < 0 ? 0 : 1);
        const Real HplusY = std::fabs(distPy)<eps ? (Real)0.5 :
                                                    (distPy < 0 ? 0 : 1);
        const Real HminuY = std::fabs(distMy)<eps ? (Real)0.5 :
                                                    (distMy < 0 ? 0 : 1);
        const Real HplusZ = std::fabs(distPz)<eps ? (Real)0.5 :
                                                    (distPz < 0 ? 0 : 1);
        const Real HminuZ = std::fabs(distMz)<eps ? (Real)0.5 :
                                                    (distMz < 0 ? 0 : 1);
        const Real gradUX = inv2h*(distPx - distMx);
        const Real gradUY = inv2h*(distPy - distMy);
        const Real gradUZ = inv2h*(distPz - distMz);
        const Real gradIX = inv2h*(IplusX - IminuX);
        const Real gradIY = inv2h*(IplusY - IminuY);
        const Real gradIZ = inv2h*(IplusZ - IminuZ);
        const Real gradHX = (HplusX - HminuX);
        const Real gradHY = (HplusY - HminuY);
        const Real gradHZ = (HplusZ - HminuZ);
        const Real gradUSq = gradUX*gradUX +gradUY*gradUY +gradUZ*gradUZ +eps;
        const Real numH    = gradIX*gradUX +gradIY*gradUY +gradIZ*gradUZ;
        const Real numD    = gradHX*gradUX +gradHY*gradUY +gradHZ*gradUZ;
        const Real Delta   = numD / gradUSq;
        const Real H       = numH / gradUSq;

        o->write(ix, iy, iz, H, Delta, gradUX, gradUY, gradUZ, fac1);
        b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
        #ifndef NDEBUG
        if(H<0||H>1+eps) printf("invalid H?: %f %f %f: %10.10e\n",x,y,z,H);
        #endif
      }
    }
  }
};
}

#if 0
namespace TorusObstacle
{
struct FillBlocks
{
  const Real big_r, small_r, smoothing_length;
  const Real position[3];
  Real sphere_box[3][2];

  void _find_sphere_box()
  {
    sphere_box[0][0] = position[0] - 2*(small_r + smoothing_length);
    sphere_box[0][1] = position[0] + 2*(small_r + smoothing_length);
    sphere_box[1][0] = position[1] - 2*(big_r + small_r + smoothing_length);
    sphere_box[1][1] = position[1] + 2*(big_r + small_r + smoothing_length);
    sphere_box[2][0] = position[2] - 2*(big_r + small_r + smoothing_length);
    sphere_box[2][1] = position[2] + 2*(big_r + small_r + smoothing_length);
  }

  FillBlocks(Real _big_r, Real _small_r, Real _smoothing_length,  Real p[3]):
    big_r(_big_r), small_r(_small_r), smoothing_length(_smoothing_length),
    position{p[0],p[1],p[2]}
  {

    _find_sphere_box();
  }

  FillBlocks(const FillBlocks& c):
    big_r(c.big_r), small_r(c.small_r), smoothing_length(c.smoothing_length),
    position{c.position[0],c.position[1],c.position[2]}
  {
    _find_sphere_box();
  }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    Real intersection[3][2] = {
        max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
        max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
        max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
  {
    Real min_pos[3], max_pos[3];

    info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i)
    {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return _is_touching(min_pos,max_pos);
  }

  Real mollified_heaviside(const Real x, const Real eps) const
  {
    const Real alpha = M_PI*min(1., max(0., (x+0.5*eps)/eps));

    return 0.5+0.5*cos(alpha);
  }


  inline void operator()(const BlockInfo& info, FluidBlock& b) const
  {
    if(_is_touching(info))
      //  if(true)
    {
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
          {
            Real p[3];
            info.pos(p, ix, iy, iz);

            const Real t[3] = {
                (Real)(p[0] - position[0]),
                (Real)(p[1] - position[1]),
                (Real)(p[2] - position[2])
            };

            //mappele //sqrt(t[1] * t[1] + t[2] * t[2]);//
            const Real r = sqrt(t[1]*t[1] + t[2]*t[2]);

            if (r > 0)
            {
              const Real c[3] = {
                  0,
                  big_r*t[1]/r,
                  big_r*t[2]/r
              };

              const Real d = sqrt(pow(t[0]-c[0], 2) +  pow(t[1]-c[1], 2) + pow(t[2]-c[2], 2));
              const Real chi = mollified_heaviside(d-small_r, smoothing_length);
              b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
            }
          }
    }
  }
};
}

namespace EllipsoidObstacle
{
// code from http://www.geometrictools.com/

static Real DistancePointEllipseSpecial (const Real e[2], const Real y[2], Real x[2])
{
  Real distance = (Real)0;
  if (y[1] > (Real)0)
  {
    if (y[0] > (Real)0)
    {
      // Bisect to compute the root of F(t) for t >= -e1*e1.
      Real esqr[2] = { e[0]*e[0], e[1]*e[1] };
      Real ey[2] = { e[0]*y[0], e[1]*y[1] };
      Real t0 = -esqr[1] + ey[1];
      Real t1 = -esqr[1] + sqrt(ey[0]*ey[0] + ey[1]*ey[1]);
      Real t = t0;
      const int imax = 2*std::numeric_limits<Real>::max_exponent;
      static constexpr Real eps = std::numeric_limits<Real>::epsilon();
      for (int i = 0; i < imax; ++i)
      {
        t = ((Real)0.5)*(t0 + t1);
        if (std::fabs(t-t0)<eps || std::fabs(t-t1)<eps)
        {
          break;
        }

        Real r[2] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]) };
        Real f = r[0]*r[0] + r[1]*r[1] - (Real)1;
        if (f > (Real)0)
        {
          t0 = t;
        }
        else if (f < (Real)0)
        {
          t1 = t;
        }
        else
        {
          break;
        }
      }

      x[0] = esqr[0]*y[0]/(t + esqr[0]);
      x[1] = esqr[1]*y[1]/(t + esqr[1]);
      Real d[2] = { x[0] - y[0], x[1] - y[1] };
      distance = sqrt(d[0]*d[0] + d[1]*d[1]);
    }
    else  // y0 == 0
    {
      x[0] = (Real)0;
      x[1] = e[1];
      distance = fabs(y[1] - e[1]);
    }
  }
  else  // y1 == 0
  {
    Real denom0 = e[0]*e[0] - e[1]*e[1];
    Real e0y0 = e[0]*y[0];
    if (e0y0 < denom0)
    {
      // y0 is inside the subinterval.
      Real x0de0 = e0y0/denom0;
      Real x0de0sqr = x0de0*x0de0;
      x[0] = e[0]*x0de0;
      x[1] = e[1]*sqrt(fabs((Real)1 - x0de0sqr));
      Real d0 = x[0] - y[0];
      distance = sqrt(d0*d0 + x[1]*x[1]);
    }
    else
    {
      // y0 is outside the subinterval.  The closest ellipse point has
      // x1 == 0 and is on the domain-boundary interval (x0/e0)^2 = 1.
      x[0] = e[0];
      x[1] = (Real)0;
      distance = fabs(y[0] - e[0]);
    }
  }
  return distance;
}
//----------------------------------------------------------------------------
// The ellipse is (x0/e0)^2 + (x1/e1)^2 = 1.  The query point is (y0,y1).
// The function returns the distance from the query point to the ellipse.
// It also computes the ellipse point (x0,x1) that is closest to (y0,y1).
//----------------------------------------------------------------------------

static Real DistancePointEllipse (const Real e[2], const Real y[2], Real x[2])
{
  // Determine reflections for y to the first quadrant.
  bool reflect[2];
  int i, j;
  for (i = 0; i < 2; ++i)
  {
    reflect[i] = (y[i] < (Real)0);
  }

  // Determine the axis order for decreasing extents.
  int permute[2];
  if (e[0] < e[1])
  {
    permute[0] = 1;  permute[1] = 0;
  }
  else
  {
    permute[0] = 0;  permute[1] = 1;
  }

  int invpermute[2];
  for (i = 0; i < 2; ++i)
  {
    invpermute[permute[i]] = i;
  }

  Real locE[2], locY[2];
  for (i = 0; i < 2; ++i)
  {
    j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if (reflect[j])
    {
      locY[i] = -locY[i];
    }
  }

  Real locX[2];
  Real distance = DistancePointEllipseSpecial(locE, locY, locX);

  // Restore the axis order and reflections.
  for (i = 0; i < 2; ++i)
  {
    j = invpermute[i];
    if (reflect[j])
    {
      locX[j] = -locX[j];
    }
    x[i] = locX[j];
  }

  return distance;
}
//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1 with e0 >= e1 >= e2.
// The query point is (y0,y1,y2) with y0 >= 0, y1 >= 0, and y2 >= 0.  The
// function returns the distance from the query point to the ellipsoid.  It
// also computes the ellipsoid point (x0,x1,x2) in the first octant that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------

static Real DistancePointEllipsoidSpecial (const Real e[3], const Real y[3],
    Real x[3])
{
  Real distance;
  if (y[2] > (Real)0)
  {
    if (y[1] > (Real)0)
    {
      if (y[0] > (Real)0)
      {
        // Bisect to compute the root of F(t) for t >= -e2*e2.
        Real esqr[3] = { e[0]*e[0], e[1]*e[1], e[2]*e[2] };
        Real ey[3] = { e[0]*y[0], e[1]*y[1], e[2]*y[2] };
        Real t0 = -esqr[2] + ey[2];
        Real t1 = -esqr[2] + sqrt(ey[0]*ey[0] + ey[1]*ey[1] +
            ey[2]*ey[2]);
        Real t = t0;
        const int imax = 2*std::numeric_limits<Real>::max_exponent;
        for (int i = 0; i < imax; ++i)
        {
          t = ((Real)0.5)*(t0 + t1);
          static constexpr Real eps = std::numeric_limits<Real>::epsilon();
          if (std::fabs(t-t0)<eps || std::fabs(t-t1)<eps)
          {
            break;
          }

          Real r[3] = { ey[0]/(t + esqr[0]), ey[1]/(t + esqr[1]),
              ey[2]/(t + esqr[2]) };
          Real f = r[0]*r[0] + r[1]*r[1] + r[2]*r[2] - (Real)1;
          if (f > (Real)0)
          {
            t0 = t;
          }
          else if (f < (Real)0)
          {
            t1 = t;
          }
          else
          {
            break;
          }
        }

        x[0] = esqr[0]*y[0]/(t + esqr[0]);
        x[1] = esqr[1]*y[1]/(t + esqr[1]);
        x[2] = esqr[2]*y[2]/(t + esqr[2]);
        Real d[3] = { x[0] - y[0], x[1] - y[1], x[2] - y[2] };
        distance = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      }
      else  // y0 == 0
      {
        x[0] = (Real)0;
        Real etmp[2] = { e[1], e[2] };
        Real ytmp[2] = { y[1], y[2] };
        Real xtmp[2];
        distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
        x[1] = xtmp[0];
        x[2] = xtmp[1];
      }
    }
    else  // y1 == 0
    {
      x[1] = (Real)0;
      if (y[0] > (Real)0)
      {
        Real etmp[2] = { e[0], e[2] };
        Real ytmp[2] = { y[0], y[2] };
        Real xtmp[2];
        distance = DistancePointEllipseSpecial(etmp, ytmp, xtmp);
        x[0] = xtmp[0];
        x[2] = xtmp[1];
      }
      else  // y0 == 0
      {
        x[0] = (Real)0;
        x[2] = e[2];
        distance = fabs(y[2] - e[2]);
      }
    }
  }
  else  // y2 == 0
  {
    Real denom[2] = { e[0]*e[0] - e[2]*e[2], e[1]*e[1] - e[2]*e[2] };
    Real ey[2] = { e[0]*y[0], e[1]*y[1] };
    if (ey[0] < denom[0] && ey[1] < denom[1])
    {
      // (y0,y1) is inside the axis-aligned bounding rectangle of the
      // subellipse.  This intermediate test is designed to guard
      // against the division by zero when e0 == e2 or e1 == e2.
      Real xde[2] = { ey[0]/denom[0], ey[1]/denom[1] };
      Real xdesqr[2] = { xde[0]*xde[0], xde[1]*xde[1] };
      Real discr = (Real)1 - xdesqr[0] - xdesqr[1];
      if (discr > (Real)0)
      {
        // (y0,y1) is inside the subellipse.  The closest ellipsoid
        // point has x2 > 0.
        x[0] = e[0]*xde[0];
        x[1] = e[1]*xde[1];
        x[2] = e[2]*sqrt(discr);
        Real d[2] = { x[0] - y[0], x[1] - y[1] };
        distance = sqrt(d[0]*d[0] + d[1]*d[1] + x[2]*x[2]);
      }
      else
      {
        // (y0,y1) is outside the subellipse.  The closest ellipsoid
        // point has x2 == 0 and is on the domain-boundary ellipse
        // (x0/e0)^2 + (x1/e1)^2 = 1.
        x[2] = (Real)0;
        distance = DistancePointEllipseSpecial(e, y, x);
      }
    }
    else
    {
      // (y0,y1) is outside the subellipse.  The closest ellipsoid
      // point has x2 == 0 and is on the domain-boundary ellipse
      // (x0/e0)^2 + (x1/e1)^2 = 1.
      x[2] = (Real)0;
      distance = DistancePointEllipseSpecial(e, y, x);
    }
  }
  return distance;
}
//----------------------------------------------------------------------------
// The ellipsoid is (x0/e0)^2 + (x1/e1)^2 + (x2/e2)^2 = 1.  The query point is
// (y0,y1,y2).  The function returns the distance from the query point to the
// ellipsoid.   It also computes the ellipsoid point (x0,x1,x2) that is
// closest to (y0,y1,y2).
//----------------------------------------------------------------------------
template <typename Real>
static Real DistancePointEllipsoid (const Real e[3], const Real y[3], Real x[3])
{
  // Determine reflections for y to the first octant.
  bool reflect[3];
  int i, j;
  for (i = 0; i < 3; ++i)
  {
    reflect[i] = (y[i] < (Real)0);
  }

  // Determine the axis order for decreasing extents.
  int permute[3];
  if (e[0] < e[1])
  {
    if (e[2] < e[0])
    {
      permute[0] = 1;  permute[1] = 0;  permute[2] = 2;
    }
    else if (e[2] < e[1])
    {
      permute[0] = 1;  permute[1] = 2;  permute[2] = 0;
    }
    else
    {
      permute[0] = 2;  permute[1] = 1;  permute[2] = 0;
    }
  }
  else
  {
    if (e[2] < e[1])
    {
      permute[0] = 0;  permute[1] = 1;  permute[2] = 2;
    }
    else if (e[2] < e[0])
    {
      permute[0] = 0;  permute[1] = 2;  permute[2] = 1;
    }
    else
    {
      permute[0] = 2;  permute[1] = 0;  permute[2] = 1;
    }
  }

  int invpermute[3];
  for (i = 0; i < 3; ++i)
  {
    invpermute[permute[i]] = i;
  }

  Real locE[3], locY[3];
  for (i = 0; i < 3; ++i)
  {
    j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if (reflect[j])
    {
      locY[i] = -locY[i];
    }
  }

  Real locX[3];
  Real distance = DistancePointEllipsoidSpecial(locE, locY, locX);

  // Restore the axis order and reflections.
  for (i = 0; i < 3; ++i)
  {
    j = invpermute[i];
    if (reflect[j])
    {
      locX[j] = -locX[j];
    }
    x[i] = locX[j];
  }

  return distance;
}

struct FillBlocksEllipsoid
{
  const Real e0,e1,e2,smoothing_length;
  const Real position[3];
  const Real quaternion[4];
  Real sphere_box[3][2];

  void _find_sphere_box()
  {
    const Real maxAxis = std::max(std::max(e0,e1),e2);
    sphere_box[0][0] = position[0] - 2*(maxAxis + smoothing_length);
    sphere_box[0][1] = position[0] + 2*(maxAxis + smoothing_length);
    sphere_box[1][0] = position[1] - 2*(maxAxis + smoothing_length);
    sphere_box[1][1] = position[1] + 2*(maxAxis + smoothing_length);
    sphere_box[2][0] = position[2] - 2*(maxAxis + smoothing_length);
    sphere_box[2][1] = position[2] + 2*(maxAxis + smoothing_length);
  }

  FillBlocksEllipsoid(const Real _e0, const Real _e1, const Real _e2, const Real sl, const Real p[3], const Real q[4]):
    e0(_e0),e1(_e1),e2(_e2), smoothing_length(sl), position{p[0],p[1],p[2]},
    quaternion{q[0],q[1],q[2],q[3]}
  {
    _find_sphere_box();
  }

  FillBlocksEllipsoid(const FillBlocksEllipsoid& c):
    e0(c.e0),e1(c.e1),e2(c.e2), smoothing_length(c.smoothing_length),
    position{c.position[0],c.position[1],c.position[2]},
    quaternion{c.quaternion[0],c.quaternion[1],c.quaternion[2],c.quaternion[3]}
  {
    _find_sphere_box();
  }

  Real mollified_heaviside(const Real x, const Real eps) const
  {
    const Real alpha = M_PI*std::min(1., std::max(0., (x+0.5*eps)/eps));

    return 0.5+0.5*cos(alpha);
  }

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    using std::min;
    using std::max;
    Real intersection[3][2] = {
        max(min_pos[0], sphere_box[0][0]), min(max_pos[0], sphere_box[0][1]),
        max(min_pos[1], sphere_box[1][0]), min(max_pos[1], sphere_box[1][1]),
        max(min_pos[2], sphere_box[2][0]), min(max_pos[2], sphere_box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;

  }

  bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
  {
    Real min_pos[3], max_pos[3];

    info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i)
    {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return _is_touching(min_pos,max_pos);
  }

  inline void operator()(const BlockInfo& info, FluidBlock& b) const
  {

    const Real w = quaternion[0];
    const Real x = quaternion[1];
    const Real y = quaternion[2];
    const Real z = quaternion[3];

    const Real Rmatrix[3][3] = {
        {1-2*(y*y+z*z),  2*(x*y+z*w),    2*(x*z-y*w)},
        {2*(x*y-z*w),    1-2*(x*x+z*z),  2*(y*z+x*w)},
        {2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x*x+y*y)}
    };

    if(_is_touching(info))
    {
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
          {
            Real p[3];
            info.pos(p, ix, iy, iz);

            // translate
            p[0] -= position[0];
            p[1] -= position[1];
            p[2] -= position[2];

            // rotate
            const Real t[3] = {
                Rmatrix[0][0]*p[0] + Rmatrix[0][1]*p[1] + Rmatrix[0][2]*p[2],
                Rmatrix[1][0]*p[0] + Rmatrix[1][1]*p[1] + Rmatrix[1][2]*p[2],
                Rmatrix[2][0]*p[0] + Rmatrix[2][1]*p[1] + Rmatrix[2][2]*p[2]
            };

            // find distance
            const Real e[3] = {e0,e1,e2};
            Real xs[3];
            const Real dist=DistancePointEllipsoid (e, t, xs);
            const int sign = ( (t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) > (xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2]) ) ? 1 : -1;
            const Real chi =  mollified_heaviside(sign*dist, smoothing_length);
            b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
          }
    }
  }
};


struct FillBlocks_Towers
{
  const Real e0,e1,e2,safe_distance;
  Real position[3];
  Real quaternion[4];
  Real normalI[3], normalJ[3], normalK[3];
  Real rotationMatrix[3][3];

  void normalizeNormals()
  {
    const Real magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
    const Real magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
    const Real magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
    assert(magI > std::numeric_limits<Real>::epsilon());
    assert(magJ > std::numeric_limits<Real>::epsilon());
    assert(magK > std::numeric_limits<Real>::epsilon());
    const Real invMagI = 1.0/magI;
    const Real invMagJ = 1.0/magJ;
    const Real invMagK = 1.0/magK;

    for(int i=0;i<3;++i)
    {
      normalI[i]=std::abs(normalI[i])*invMagI;
      normalJ[i]=std::abs(normalJ[i])*invMagJ;
      normalK[i]=std::abs(normalK[i])*invMagK;
    }
  }

  void _computeRotationMatrix()
  {
    const Real w = quaternion[0];
    const Real x = quaternion[1];
    const Real y = quaternion[2];
    const Real z = quaternion[3];

    const Real Rmatrix[3][3] = {
        {1-2*(y*y+z*z),  2*(x*y+z*w),    2*(x*z-y*w)},
        {2*(x*y-z*w),    1-2*(x*x+z*z),  2*(y*z+x*w)},
        {2*(x*z+y*w),    2*(y*z-x*w),    1-2*(x*x+y*y)}
    };
    memcpy(rotationMatrix, Rmatrix, sizeof(Rmatrix));

    // initial axes are grid-aligned: e0 along x, e1 along y, e2 along z
    const Real nx[3] = {1.0,0.0,0.0};
    const Real ny[3] = {0.0,1.0,0.0};
    const Real nz[3] = {0.0,0.0,1.0};

    // now we rotate to computational frame
    for(int i=0;i<3;++i)
    {
      normalI[i] = Rmatrix[i][0]*nx[0] + Rmatrix[i][1]*nx[1] + Rmatrix[i][2]*nx[2];
      normalJ[i] = Rmatrix[i][0]*ny[0] + Rmatrix[i][1]*ny[1] + Rmatrix[i][2]*ny[2];
      normalK[i] = Rmatrix[i][0]*nz[0] + Rmatrix[i][1]*nz[1] + Rmatrix[i][2]*nz[2];
    }

    normalizeNormals();
  }



  FillBlocks_Towers(const Real _e0, const Real _e1, const Real _e2, const Real max_dx, const Real p[3], const Real q[4]):
    e0(_e0),e1(_e1),e2(_e2), safe_distance(max_dx), position{p[0],p[1],p[2]},
    quaternion{q[0],q[1],q[2],q[3]}
  {
    _computeRotationMatrix();
  }
  FillBlocks_Towers(const FillBlocks_Towers& c):
    e0(c.e0),e1(c.e1),e2(c.e2), safe_distance(c.safe_distance),
    position{c.position[0],c.position[1],c.position[2]},
    quaternion{c.quaternion[0],c.quaternion[1],c.quaternion[2],c.quaternion[3]}
  {
    _computeRotationMatrix();
  }

  bool _is_touching(const Real start[3], const Real end[3]) const
  {
    const Real AABB_w[3] = {
        (Real)0.5*(end[0] - start[0]),
        (Real)0.5*(end[1] - start[1]),
        (Real)0.5*(end[2] - start[2])
    }; // halfwidth

    const Real AABB_c[3] = {
        start[0] + AABB_w[0],
        start[1] + AABB_w[1],
        start[2] + AABB_w[2]
    }; // center

    assert(AABB_w[0]>0);
    assert(AABB_w[1]>0);
    assert(AABB_w[2]>0);

    const Real w[3] = {e0,e1,e2};
    assert(w[0]>0);
    assert(w[1]>0);
    assert(w[2]>0);

    bool intersects = true;
    Real r;
    {
      r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
      intersects &= ((position[0]-r <= AABB_c[0] + AABB_w[0]) && (position[0]+r >= AABB_c[0] - AABB_w[0]));

      r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
      intersects &= ((position[1]-r <= AABB_c[1] + AABB_w[1]) && (position[1]+r >= AABB_c[1] - AABB_w[1]));

      r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
      intersects &= ((position[2]-r <= AABB_c[2] + AABB_w[2]) && (position[2]+r >= AABB_c[2] - AABB_w[2]));
    }
    {
      r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
      intersects &= ((AABB_c[0]-r <= position[0] + w[0]) && (AABB_c[0]+r >= position[0] - w[0]));

      r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
      intersects &= ((AABB_c[1]-r <= position[1] + w[1]) && (AABB_c[1]+r >= position[1] - w[1]));

      r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
      intersects &= ((AABB_c[2]-r <= position[2] + w[2]) && (AABB_c[2]+r >= position[2] - w[2]));
    }
    return intersects;
  }

  bool _is_touching(const BlockInfo& info, const int buffer_dx=0) const
  {

    const Real safe_distance_info = 2.0*info.h_gridpoint; // two points on each side needed for towers

    // AABB - OBB intersection.
    // ellipsoid is inside OBB with normals normalI,normalJ,normalK and widths e0,e1,e2
    // grid block is AABB

    Real start[3], end[3];

    info.pos(start, 0,0,0);
    info.pos(end, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);

    for(int i=0;i<3;++i)
    {
      start[i]-=(safe_distance_info + buffer_dx*info.h_gridpoint);
      end[i] += (safe_distance_info + buffer_dx*info.h_gridpoint);
    }

    return _is_touching(start,end);
  }


  inline Real distanceToEllipsoid(const Real x, const Real y, const Real z) const
  {
    const Real e[3] = {e0,e1,e2};
    const Real t[3] = {x,y,z};
    Real xs[3];
    const Real dist=DistancePointEllipsoid (e, t, xs);
    const int invSign = ( (t[0]*t[0]+t[1]*t[1]+t[2]*t[2]) > (xs[0]*xs[0]+xs[1]*xs[1]+xs[2]*xs[2]) ) ? -1 : 1;
    return dist*invSign;// pos inside, neg outside
  }

  Real getHeavisideFDMH1(const Real x, const Real y, const Real z, const Real h) const
  {
    const Real dist = distanceToEllipsoid(x,y,z);
    if(dist >= +h) return 1;
    if(dist <= -h) return 0;
    assert(std::abs(dist)<=h);

    const Real distPx = distanceToEllipsoid(x+h,y,z);
    const Real distMx = distanceToEllipsoid(x-h,y,z);
    const Real distPy = distanceToEllipsoid(x,y+h,z);
    const Real distMy = distanceToEllipsoid(x,y-h,z);
    const Real distPz = distanceToEllipsoid(x,y,z+h);
    const Real distMz = distanceToEllipsoid(x,y,z-h);

    // compute first primitive of H(x): I(x) = int_0^x H(y) dy and set it to zero outside the cylinder
    const Real IplusX = distPx < 0 ? 0 : distPx;
    const Real IminuX = distMx < 0 ? 0 : distMx;
    const Real IplusY = distPy < 0 ? 0 : distPy;
    const Real IminuY = distMy < 0 ? 0 : distMy;
    const Real IplusZ = distPz < 0 ? 0 : distPz;
    const Real IminuZ = distMz < 0 ? 0 : distMz;

    assert(IplusX>=0);
    assert(IminuX>=0);
    assert(IplusY>=0);
    assert(IminuY>=0);
    assert(IplusZ>=0);
    assert(IminuZ>=0);


    // gradI
    const Real gradIX = 0.5/h * (IplusX - IminuX);
    const Real gradIY = 0.5/h * (IplusY - IminuY);
    const Real gradIZ = 0.5/h * (IplusZ - IminuZ);

    // gradU
    const Real gradUX = 0.5/h * (distPx - distMx);
    const Real gradUY = 0.5/h * (distPy - distMy);
    const Real gradUZ = 0.5/h * (distPz - distMz);

    const Real H = (gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ)/(gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ);

    //        assert(H>=0 && H<=1);
    if(H<0.0 || H>1.0)
      printf("invalid H?: %10.10e %10.10e %10.10e: %10.10e\n",x,y,z,H);

    return H;
  }


  inline void operator()(const BlockInfo& info, FluidBlock& b) const
  {
    if(_is_touching(info))
    {
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
        for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
          {
            Real p[3];
            info.pos(p, ix, iy, iz);

            // translate
            p[0] -= position[0];
            p[1] -= position[1];
            p[2] -= position[2];

            // rotate
            const Real t[3] = {
                rotationMatrix[0][0]*p[0] + rotationMatrix[0][1]*p[1] + rotationMatrix[0][2]*p[2],
                rotationMatrix[1][0]*p[0] + rotationMatrix[1][1]*p[1] + rotationMatrix[1][2]*p[2],
                rotationMatrix[2][0]*p[0] + rotationMatrix[2][1]*p[1] + rotationMatrix[2][2]*p[2]
            };

            const Real chi = getHeavisideFDMH1(t[0],t[1],t[2],info.h_gridpoint);
            b(ix, iy, iz).chi = std::max(chi, b(ix, iy, iz).chi);
          }
    }
  }
};

}
#endif

#endif
