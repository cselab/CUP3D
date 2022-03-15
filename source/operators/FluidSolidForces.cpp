//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "FluidSolidForces.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

struct KernelComputeForces : public ObstacleVisitor
{
  ObstacleVector * const obstacle_vector;
  const Real nu, dt;
  LabMPI * lab_ptr = nullptr;
  const BlockInfo * info_ptr = nullptr;

  const int big   = 5;
  const int small = -4;

  const int stencil_start[3] = {small,small,small}, stencil_end[3] = {big,big,big};
  StencilInfo stencil{small,small,small, big,big,big, true, {FE_CHI,FE_U,FE_V,FE_W}};

  KernelComputeForces(Real _nu, Real _dt, ObstacleVector* ov) :
    obstacle_vector(ov), nu(_nu), dt(_dt) { }

  void operator()(LabMPI&lab,const BlockInfo&info)
  {
    // first store the lab and info, then do visitor
    lab_ptr = & lab;
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    lab_ptr = nullptr;
    info_ptr = nullptr;
  }

  void visit(Obstacle* const op)
  {

    LabMPI& l = * lab_ptr;
    const BlockInfo& info = * info_ptr;
    assert(lab_ptr not_eq nullptr && info_ptr not_eq nullptr);
    const std::vector<ObstacleBlock*>& obstblocks = op->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;
    if (o->nPoints == 0) return;
    assert(o->filled);

    const std::array<Real,3> CM = op->getCenterOfMass();
    const std::array<Real,3> omega = op->getAngularVelocity();
    const std::array<Real,3> uTrans = op->getTranslationVelocity();
    Real velUnit[3] = {0., 0., 0.};
    const Real vel_norm = std::sqrt(uTrans[0]*uTrans[0]
                                    + uTrans[1]*uTrans[1]
                                    + uTrans[2]*uTrans[2]);
    if (vel_norm>1e-9) {
        velUnit[0] = uTrans[0] / vel_norm;
        velUnit[1] = uTrans[1] / vel_norm;
        velUnit[2] = uTrans[2] / vel_norm;
    }

    const Real _1oH = nu / info.h;
    //const Real _h3 = std::pow(info.h,3);

    //loop over elements of block info that have nonzero gradChi
    for(int i=0; i<o->nPoints; i++)
    {
      Real p[3];
      const int ix = o->surface[i]->ix;
      const int iy = o->surface[i]->iy;
      const int iz = o->surface[i]->iz;

      // WAS A SOURCE OF HUGE BUG - due to unzeroed values in surfData arrays. The kid had forgotten to initialize allocated arrays in surfData to zero!!
      //if(o->chi[iz][iy][ix] < 1e-16) continue;

      info.pos(p, ix, iy, iz);

      //shear stresses
      const Real normX = o->surface[i]->dchidx; //*h^3 (multiplied in dchidx)
      const Real normY = o->surface[i]->dchidy; //*h^3 (multiplied in dchidy)
      const Real normZ = o->surface[i]->dchidz; //*h^3 (multiplied in dchidz)
      const Real norm = 1.0/std::sqrt(normX*normX+normY*normY+normZ*normZ);

      Real dx = normX*norm;
      Real dy = normY*norm;
      Real dz = normZ*norm;

      int x = ix;
      int y = iy;
      int z = iz;
      int found = 0;
      for (int kk = 1 ; kk < 10 ; kk++) //10 is arbitrary
      {
        if ((int)abs(kk*dx) > 3 || (int)abs(kk*dy) > 3 || (int)abs(kk*dz) > 3) break; //3 means we moved too far
        if (l(x,y,z).chi <0.3 && found == 0) break;
        x  = ix + kk*dx; 
        y  = iy + kk*dy;
        z  = iz + kk*dz;
        if (l(x,y,z).chi < 0.3 ) found ++;
      }

      Real dudx1 = normX > 0 ? (l(x+1,y,z).u-l(x,y,z).u) : (l(x,y,z).u-l(x-1,y,z).u);
      Real dvdx1 = normX > 0 ? (l(x+1,y,z).v-l(x,y,z).v) : (l(x,y,z).v-l(x-1,y,z).v);
      Real dwdx1 = normX > 0 ? (l(x+1,y,z).w-l(x,y,z).w) : (l(x,y,z).w-l(x-1,y,z).w);
      Real dudy1 = normY > 0 ? (l(x,y+1,z).u-l(x,y,z).u) : (l(x,y,z).u-l(x,y-1,z).u);
      Real dvdy1 = normY > 0 ? (l(x,y+1,z).v-l(x,y,z).v) : (l(x,y,z).v-l(x,y-1,z).v);
      Real dwdy1 = normY > 0 ? (l(x,y+1,z).w-l(x,y,z).w) : (l(x,y,z).w-l(x,y-1,z).w);
      Real dudz1 = normY > 0 ? (l(x,y,z+1).u-l(x,y,z).u) : (l(x,y,z).u-l(x,y,z-1).u);
      Real dvdz1 = normY > 0 ? (l(x,y,z+1).v-l(x,y,z).v) : (l(x,y,z).v-l(x,y,z-1).v);
      Real dwdz1 = normY > 0 ? (l(x,y,z+1).w-l(x,y,z).w) : (l(x,y,z).w-l(x,y,z-1).w);
      Real dudx2 = 0.0;
      Real dvdx2 = 0.0;
      Real dwdx2 = 0.0;
      Real dudy2 = 0.0;
      Real dvdy2 = 0.0;
      Real dwdy2 = 0.0;
      Real dudz2 = 0.0;
      Real dvdz2 = 0.0;
      Real dwdz2 = 0.0;

      Real dudx3 = 0.0;
      Real dvdx3 = 0.0;
      Real dwdx3 = 0.0;
      Real dudy3 = 0.0;
      Real dvdy3 = 0.0;
      Real dwdy3 = 0.0;
      Real dudz3 = 0.0;
      Real dvdz3 = 0.0;
      Real dwdz3 = 0.0;

      Real dudxdy1 = 0.0;
      Real dvdxdy1 = 0.0;
      Real dwdxdy1 = 0.0;
      Real dudxdz1 = 0.0;
      Real dvdxdz1 = 0.0;
      Real dwdxdz1 = 0.0;
      Real dudydz1 = 0.0;
      Real dvdydz1 = 0.0;
      Real dwdydz1 = 0.0;

      dudxdy1 = 0.25*(l(x+1,y+1,z).u+l(x-1,y-1,z).u-l(x+1,y-1,z).u-l(x-1,y+1,z).u);
      dvdxdy1 = 0.25*(l(x+1,y+1,z).v+l(x-1,y-1,z).v-l(x+1,y-1,z).v-l(x-1,y+1,z).v);
      dwdxdy1 = 0.25*(l(x+1,y+1,z).w+l(x-1,y-1,z).w-l(x+1,y-1,z).w-l(x-1,y+1,z).w);
      if (normX > 0 && normY > 0)
      {
        dudxdy1 = (l(x+1,y+1,z).u+l(x,y,z).u-l(x+1,y,z).u-l(x,y+1,z).u);
        dvdxdy1 = (l(x+1,y+1,z).v+l(x,y,z).v-l(x+1,y,z).v-l(x,y+1,z).v);
        dwdxdy1 = (l(x+1,y+1,z).w+l(x,y,z).w-l(x+1,y,z).w-l(x,y+1,z).w);
        if ((x+2 < big) && (y+2 < big))
        {
           dudxdy1 = -0.5*( -1.5*l(x+2,y,z).u+2*l(x+2,y+1,z).u-0.5*l(x+2,y+2,z).u ) + 2*(-1.5*l(x+1,y,z).u+2*l(x+1,y+1,z).u-0.5*l(x+1,y+2,z).u) -1.5*(-1.5*l(x,y,z).u+2*l(x,y+1,z).u-0.5*l(x,y+2,z).u);
           dvdxdy1 = -0.5*( -1.5*l(x+2,y,z).v+2*l(x+2,y+1,z).v-0.5*l(x+2,y+2,z).v ) + 2*(-1.5*l(x+1,y,z).v+2*l(x+1,y+1,z).v-0.5*l(x+1,y+2,z).v) -1.5*(-1.5*l(x,y,z).v+2*l(x,y+1,z).v-0.5*l(x,y+2,z).v);
           dwdxdy1 = -0.5*( -1.5*l(x+2,y,z).w+2*l(x+2,y+1,z).w-0.5*l(x+2,y+2,z).w ) + 2*(-1.5*l(x+1,y,z).w+2*l(x+1,y+1,z).w-0.5*l(x+1,y+2,z).w) -1.5*(-1.5*l(x,y,z).w+2*l(x,y+1,z).w-0.5*l(x,y+2,z).w);
        }
      }
      if (normX < 0 && normY > 0)
      {
        dudxdy1 = (l(x,y+1,z).u+l(x-1,y,z).u-l(x,y,z).u-l(x-1,y+1,z).u);
        dvdxdy1 = (l(x,y+1,z).v+l(x-1,y,z).v-l(x,y,z).v-l(x-1,y+1,z).v);
        dwdxdy1 = (l(x,y+1,z).w+l(x-1,y,z).w-l(x,y,z).w-l(x-1,y+1,z).w);
        if ((y+2 < big) && (x-2 >= small))
        {
           dudxdy1 = 0.5*( -1.5*l(x-2,y,z).u+2*l(x-2,y+1,z).u-0.5*l(x-2,y+2,z).u ) - 2*(-1.5*l(x-1,y,z).u+2*l(x-1,y+1,z).u-0.5*l(x-1,y+2,z).u)+1.5*(-1.5*l(x,y,z).u+2*l(x,y+1,z).u-0.5*l(x,y+2,z).u);
           dvdxdy1 = 0.5*( -1.5*l(x-2,y,z).v+2*l(x-2,y+1,z).v-0.5*l(x-2,y+2,z).v ) - 2*(-1.5*l(x-1,y,z).v+2*l(x-1,y+1,z).v-0.5*l(x-1,y+2,z).v)+1.5*(-1.5*l(x,y,z).v+2*l(x,y+1,z).v-0.5*l(x,y+2,z).v);
           dwdxdy1 = 0.5*( -1.5*l(x-2,y,z).w+2*l(x-2,y+1,z).w-0.5*l(x-2,y+2,z).w ) - 2*(-1.5*l(x-1,y,z).w+2*l(x-1,y+1,z).w-0.5*l(x-1,y+2,z).w)+1.5*(-1.5*l(x,y,z).w+2*l(x,y+1,z).w-0.5*l(x,y+2,z).w);
        }
      }
      if (normX > 0 && normY < 0)
      {
        dudxdy1 = (l(x+1,y,z).u+l(x,y-1,z).u-l(x+1,y-1,z).u-l(x,y,z).u);
        dvdxdy1 = (l(x+1,y,z).v+l(x,y-1,z).v-l(x+1,y-1,z).v-l(x,y,z).v);
        dwdxdy1 = (l(x+1,y,z).w+l(x,y-1,z).w-l(x+1,y-1,z).w-l(x,y,z).w);
        if ((x+2 < big) && (y-2 >= small))
        {
           dudxdy1 = -0.5*( 1.5*l(x+2,y,z).u-2*l(x+2,y-1,z).u+0.5*l(x+2,y-2,z).u ) + 2*(1.5*l(x+1,y,z).u-2*l(x+1,y-1,z).u+0.5*l(x+1,y-2,z).u) -1.5*(1.5*l(x,y,z).u-2*l(x,y-1,z).u+0.5*l(x,y-2,z).u);
           dvdxdy1 = -0.5*( 1.5*l(x+2,y,z).v-2*l(x+2,y-1,z).v+0.5*l(x+2,y-2,z).v ) + 2*(1.5*l(x+1,y,z).v-2*l(x+1,y-1,z).v+0.5*l(x+1,y-2,z).v) -1.5*(1.5*l(x,y,z).v-2*l(x,y-1,z).v+0.5*l(x,y-2,z).v);
           dwdxdy1 = -0.5*( 1.5*l(x+2,y,z).w-2*l(x+2,y-1,z).w+0.5*l(x+2,y-2,z).w ) + 2*(1.5*l(x+1,y,z).w-2*l(x+1,y-1,z).w+0.5*l(x+1,y-2,z).w) -1.5*(1.5*l(x,y,z).w-2*l(x,y-1,z).w+0.5*l(x,y-2,z).w);
        }
      }
      if (normX < 0 && normY < 0)
      {
        dudxdy1 = (l(x,y,z).u+l(x-1,y-1,z).u-l(x,y-1,z).u-l(x-1,y,z).u);
        dvdxdy1 = (l(x,y,z).v+l(x-1,y-1,z).v-l(x,y-1,z).v-l(x-1,y,z).v);
        dwdxdy1 = (l(x,y,z).w+l(x-1,y-1,z).w-l(x,y-1,z).w-l(x-1,y,z).w);
        if ((x-2 >= small) && (y-2 >= small))
        {
           dudxdy1 = 0.5*( 1.5*l(x-2,y,z).u-2*l(x-2,y-1,z).u+0.5*l(x-2,y-2,z).u ) - 2*(1.5*l(x-1,y,z).u-2*l(x-1,y-1,z).u+0.5*l(x-1,y-2,z).u) +1.5*(1.5*l(x,y,z).u-2*l(x,y-1,z).u+0.5*l(x,y-2,z).u);
           dvdxdy1 = 0.5*( 1.5*l(x-2,y,z).v-2*l(x-2,y-1,z).v+0.5*l(x-2,y-2,z).v ) - 2*(1.5*l(x-1,y,z).v-2*l(x-1,y-1,z).v+0.5*l(x-1,y-2,z).v) +1.5*(1.5*l(x,y,z).v-2*l(x,y-1,z).v+0.5*l(x,y-2,z).v);
           dwdxdy1 = 0.5*( 1.5*l(x-2,y,z).w-2*l(x-2,y-1,z).w+0.5*l(x-2,y-2,z).w ) - 2*(1.5*l(x-1,y,z).w-2*l(x-1,y-1,z).w+0.5*l(x-1,y-2,z).w) +1.5*(1.5*l(x,y,z).w-2*l(x,y-1,z).w+0.5*l(x,y-2,z).w);
        }
      }

      dudxdz1 = 0.25*(l(x+1,y,z+1).u+l(x-1,y,z-1).u-l(x+1,y,z-1).u-l(x-1,y,z+1).u);
      dvdxdz1 = 0.25*(l(x+1,y,z+1).v+l(x-1,y,z-1).v-l(x+1,y,z-1).v-l(x-1,y,z+1).v);
      dwdxdz1 = 0.25*(l(x+1,y,z+1).w+l(x-1,y,z-1).w-l(x+1,y,z-1).w-l(x-1,y,z+1).w);
      if (normX > 0 && normZ > 0)
      {
        dudxdz1 = (l(x+1,y,z+1).u+l(x,y,z).u-l(x+1,y,z).u-l(x,y,z+1).u);
        dvdxdz1 = (l(x+1,y,z+1).v+l(x,y,z).v-l(x+1,y,z).v-l(x,y,z+1).v);
        dwdxdz1 = (l(x+1,y,z+1).w+l(x,y,z).w-l(x+1,y,z).w-l(x,y,z+1).w);
        if ((x+2 < big) && (z+2 < big))
        {
           dudxdz1 = -0.5*( -1.5*l(x+2,y,z).u+2*l(x+2,y,z+1).u-0.5*l(x+2,y,z+2).u ) + 2*(-1.5*l(x+1,y,z).u+2*l(x+1,y,z+1).u-0.5*l(x+1,y,z+2).u) -1.5*(-1.5*l(x,y,z).u+2*l(x,y,z+1).u-0.5*l(x,y,z+2).u);
           dvdxdz1 = -0.5*( -1.5*l(x+2,y,z).v+2*l(x+2,y,z+1).v-0.5*l(x+2,y,z+2).v ) + 2*(-1.5*l(x+1,y,z).v+2*l(x+1,y,z+1).v-0.5*l(x+1,y,z+2).v) -1.5*(-1.5*l(x,y,z).v+2*l(x,y,z+1).v-0.5*l(x,y,z+2).v);
           dwdxdz1 = -0.5*( -1.5*l(x+2,y,z).w+2*l(x+2,y,z+1).w-0.5*l(x+2,y,z+2).w ) + 2*(-1.5*l(x+1,y,z).w+2*l(x+1,y,z+1).w-0.5*l(x+1,y,z+2).w) -1.5*(-1.5*l(x,y,z).w+2*l(x,y,z+1).w-0.5*l(x,y,z+2).w);
        }
      }
      if (normX < 0 && normZ > 0)
      {
        dudxdz1 = (l(x,y,z+1).u+l(x-1,y,z).u-l(x,y,z).u-l(x-1,y,z+1).u);
        dvdxdz1 = (l(x,y,z+1).v+l(x-1,y,z).v-l(x,y,z).v-l(x-1,y,z+1).v);
        dwdxdz1 = (l(x,y,z+1).w+l(x-1,y,z).w-l(x,y,z).w-l(x-1,y,z+1).w);
        if ((z+2 < big) && (x-2 >= small))
        {
           dudxdz1 = 0.5*( -1.5*l(x-2,y,z).u+2*l(x-2,y,z+1).u-0.5*l(x-2,y,z+2).u ) - 2*(-1.5*l(x-1,y,z).u+2*l(x-1,y,z+1).u-0.5*l(x-1,y,z+2).u)+1.5*(-1.5*l(x,y,z).u+2*l(x,y,z+1).u-0.5*l(x,y,z+2).u);
           dvdxdz1 = 0.5*( -1.5*l(x-2,y,z).v+2*l(x-2,y,z+1).v-0.5*l(x-2,y,z+2).v ) - 2*(-1.5*l(x-1,y,z).v+2*l(x-1,y,z+1).v-0.5*l(x-1,y,z+2).v)+1.5*(-1.5*l(x,y,z).v+2*l(x,y,z+1).v-0.5*l(x,y,z+2).v);
           dwdxdz1 = 0.5*( -1.5*l(x-2,y,z).w+2*l(x-2,y,z+1).w-0.5*l(x-2,y,z+2).w ) - 2*(-1.5*l(x-1,y,z).w+2*l(x-1,y,z+1).w-0.5*l(x-1,y,z+2).w)+1.5*(-1.5*l(x,y,z).w+2*l(x,y,z+1).w-0.5*l(x,y,z+2).w);
        }
      }
      if (normX > 0 && normZ < 0)
      {
        dudxdz1 = (l(x+1,y,z).u+l(x,y,z-1).u-l(x+1,y,z-1).u-l(x,y,z).u);
        dvdxdz1 = (l(x+1,y,z).v+l(x,y,z-1).v-l(x+1,y,z-1).v-l(x,y,z).v);
        dwdxdz1 = (l(x+1,y,z).w+l(x,y,z-1).w-l(x+1,y,z-1).w-l(x,y,z).w);
        if ((x+2 < big) && (z-2 >= small))
        {
           dudxdz1 = -0.5*( 1.5*l(x+2,y,z).u-2*l(x+2,y,z-1).u+0.5*l(x+2,y,z-2).u ) + 2*(1.5*l(x+1,y,z).u-2*l(x+1,y,z-1).u+0.5*l(x+1,y,z-2).u) -1.5*(1.5*l(x,y,z).u-2*l(x,y,z-1).u+0.5*l(x,y,z-2).u);
           dvdxdz1 = -0.5*( 1.5*l(x+2,y,z).v-2*l(x+2,y,z-1).v+0.5*l(x+2,y,z-2).v ) + 2*(1.5*l(x+1,y,z).v-2*l(x+1,y,z-1).v+0.5*l(x+1,y,z-2).v) -1.5*(1.5*l(x,y,z).v-2*l(x,y,z-1).v+0.5*l(x,y,z-2).v);
           dwdxdz1 = -0.5*( 1.5*l(x+2,y,z).w-2*l(x+2,y,z-1).w+0.5*l(x+2,y,z-2).w ) + 2*(1.5*l(x+1,y,z).w-2*l(x+1,y,z-1).w+0.5*l(x+1,y,z-2).w) -1.5*(1.5*l(x,y,z).w-2*l(x,y,z-1).w+0.5*l(x,y,z-2).w);
        }
      }
      if (normX < 0 && normZ < 0)
      {
        dudxdz1 = (l(x,y,z).u+l(x-1,y,z-1).u-l(x,y,z-1).u-l(x-1,y,z).u);
        dvdxdz1 = (l(x,y,z).v+l(x-1,y,z-1).v-l(x,y,z-1).v-l(x-1,y,z).v);
        dwdxdz1 = (l(x,y,z).w+l(x-1,y,z-1).w-l(x,y,z-1).w-l(x-1,y,z).w);
        if ((x-2 >= small) && (z-2 >= small))
        {
           dudxdz1 = 0.5*( 1.5*l(x-2,y,z).u-2*l(x-2,y,z-1).u+0.5*l(x-2,y,z-2).u ) - 2*(1.5*l(x-1,y,z).u-2*l(x-1,y,z-1).u+0.5*l(x-1,y,z-2).u) +1.5*(1.5*l(x,y,z).u-2*l(x,y,z-1).u+0.5*l(x,y,z-2).u);
           dvdxdz1 = 0.5*( 1.5*l(x-2,y,z).v-2*l(x-2,y,z-1).v+0.5*l(x-2,y,z-2).v ) - 2*(1.5*l(x-1,y,z).v-2*l(x-1,y,z-1).v+0.5*l(x-1,y,z-2).v) +1.5*(1.5*l(x,y,z).v-2*l(x,y,z-1).v+0.5*l(x,y,z-2).v);
           dwdxdz1 = 0.5*( 1.5*l(x-2,y,z).w-2*l(x-2,y,z-1).w+0.5*l(x-2,y,z-2).w ) - 2*(1.5*l(x-1,y,z).w-2*l(x-1,y,z-1).w+0.5*l(x-1,y,z-2).w) +1.5*(1.5*l(x,y,z).w-2*l(x,y,z-1).w+0.5*l(x,y,z-2).w);
        }
      }

      dudydz1 = 0.25*(l(x,y+1,z+1).u+l(x,y-1,z-1).u-l(x,y+1,z-1).u-l(x,y-1,z+1).u);
      dvdydz1 = 0.25*(l(x,y+1,z+1).v+l(x,y-1,z-1).v-l(x,y+1,z-1).v-l(x,y-1,z+1).v);
      dwdydz1 = 0.25*(l(x,y+1,z+1).w+l(x,y-1,z-1).w-l(x,y+1,z-1).w-l(x,y-1,z+1).w);
      if (normY > 0 && normZ > 0)
      {
        dudydz1 = (l(x,y+1,z+1).u+l(x,y,z).u-l(x,y+1,z).u-l(x,y,z+1).u);
        dvdydz1 = (l(x,y+1,z+1).v+l(x,y,z).v-l(x,y+1,z).v-l(x,y,z+1).v);
        dwdydz1 = (l(x,y+1,z+1).w+l(x,y,z).w-l(x,y+1,z).w-l(x,y,z+1).w);
        if ((y+2 < big) && (z+2 < big))
        {
           dudydz1 = -0.5*( -1.5*l(x,y+2,z).u+2*l(x,y+2,z+1).u-0.5*l(x,y+2,z+2).u ) + 2*(-1.5*l(x,y+1,z).u+2*l(x,y+1,z+1).u-0.5*l(x,y+1,z+2).u) -1.5*(-1.5*l(x,y,z).u+2*l(x,y,z+1).u-0.5*l(x,y,z+2).u);
           dvdydz1 = -0.5*( -1.5*l(x,y+2,z).v+2*l(x,y+2,z+1).v-0.5*l(x,y+2,z+2).v ) + 2*(-1.5*l(x,y+1,z).v+2*l(x,y+1,z+1).v-0.5*l(x,y+1,z+2).v) -1.5*(-1.5*l(x,y,z).v+2*l(x,y,z+1).v-0.5*l(x,y,z+2).v);
           dwdydz1 = -0.5*( -1.5*l(x,y+2,z).w+2*l(x,y+2,z+1).w-0.5*l(x,y+2,z+2).w ) + 2*(-1.5*l(x,y+1,z).w+2*l(x,y+1,z+1).w-0.5*l(x,y+1,z+2).w) -1.5*(-1.5*l(x,y,z).w+2*l(x,y,z+1).w-0.5*l(x,y,z+2).w);
        }
      }
      if (normY < 0 && normZ > 0)
      {
        dudydz1 = (l(x,y,z+1).u+l(x,y-1,z).u-l(x,y,z).u-l(x,y-1,z+1).u);
        dvdydz1 = (l(x,y,z+1).v+l(x,y-1,z).v-l(x,y,z).v-l(x,y-1,z+1).v);
        dwdydz1 = (l(x,y,z+1).w+l(x,y-1,z).w-l(x,y,z).w-l(x,y-1,z+1).w);
        if ((z+2 < big) && (y-2 >= small))
        {
           dudydz1 = 0.5*( -1.5*l(x,y-2,z).u+2*l(x,y-2,z+1).u-0.5*l(x,y-2,z+2).u ) - 2*(-1.5*l(x,y-1,z).u+2*l(x,y-1,z+1).u-0.5*l(x,y-1,z+2).u)+1.5*(-1.5*l(x,y,z).u+2*l(x,y,z+1).u-0.5*l(x,y,z+2).u);
           dvdydz1 = 0.5*( -1.5*l(x,y-2,z).v+2*l(x,y-2,z+1).v-0.5*l(x,y-2,z+2).v ) - 2*(-1.5*l(x,y-1,z).v+2*l(x,y-1,z+1).v-0.5*l(x,y-1,z+2).v)+1.5*(-1.5*l(x,y,z).v+2*l(x,y,z+1).v-0.5*l(x,y,z+2).v);
           dwdydz1 = 0.5*( -1.5*l(x,y-2,z).w+2*l(x,y-2,z+1).w-0.5*l(x,y-2,z+2).w ) - 2*(-1.5*l(x,y-1,z).w+2*l(x,y-1,z+1).w-0.5*l(x,y-1,z+2).w)+1.5*(-1.5*l(x,y,z).w+2*l(x,y,z+1).w-0.5*l(x,y,z+2).w);
        }
      }
      if (normY > 0 && normZ < 0)
      {
        dudydz1 = (l(x,y+1,z).u+l(x,y,z-1).u-l(x,y+1,z-1).u-l(x,y,z).u);
        dvdydz1 = (l(x,y+1,z).v+l(x,y,z-1).v-l(x,y+1,z-1).v-l(x,y,z).v);
        dwdydz1 = (l(x,y+1,z).w+l(x,y,z-1).w-l(x,y+1,z-1).w-l(x,y,z).w);
        if ((y+2 < big) && (z-2 >= small))
        {
           dudydz1 = -0.5*( 1.5*l(x,y+2,z).u-2*l(x,y+2,z-1).u+0.5*l(x,y+2,z-2).u ) + 2*(1.5*l(x,y+1,z).u-2*l(x,y+1,z-1).u+0.5*l(x,y+1,z-2).u) -1.5*(1.5*l(x,y,z).u-2*l(x,y,z-1).u+0.5*l(x,y,z-2).u);
           dvdydz1 = -0.5*( 1.5*l(x,y+2,z).v-2*l(x,y+2,z-1).v+0.5*l(x,y+2,z-2).v ) + 2*(1.5*l(x,y+1,z).v-2*l(x,y+1,z-1).v+0.5*l(x,y+1,z-2).v) -1.5*(1.5*l(x,y,z).v-2*l(x,y,z-1).v+0.5*l(x,y,z-2).v);
           dwdydz1 = -0.5*( 1.5*l(x,y+2,z).w-2*l(x,y+2,z-1).w+0.5*l(x,y+2,z-2).w ) + 2*(1.5*l(x,y+1,z).w-2*l(x,y+1,z-1).w+0.5*l(x,y+1,z-2).w) -1.5*(1.5*l(x,y,z).w-2*l(x,y,z-1).w+0.5*l(x,y,z-2).w);
        }
      }
      if (normY < 0 && normZ < 0)
      {
        dudydz1 = (l(x,y,z).u+l(x,y-1,z-1).u-l(x,y,z-1).u-l(x,y-1,z).u);
        dvdydz1 = (l(x,y,z).v+l(x,y-1,z-1).v-l(x,y,z-1).v-l(x,y-1,z).v);
        dwdydz1 = (l(x,y,z).w+l(x,y-1,z-1).w-l(x,y,z-1).w-l(x,y-1,z).w);
        if ((y-2 >= small) && (z-2 >= small))
        {
           dudydz1 = 0.5*( 1.5*l(x,y-2,z).u-2*l(x,y-2,z-1).u+0.5*l(x,y-2,z-2).u ) - 2*(1.5*l(x,y-1,z).u-2*l(x,y-1,z-1).u+0.5*l(x,y-1,z-2).u) +1.5*(1.5*l(x,y,z).u-2*l(x,y,z-1).u+0.5*l(x,y,z-2).u);
           dvdydz1 = 0.5*( 1.5*l(x,y-2,z).v-2*l(x,y-2,z-1).v+0.5*l(x,y-2,z-2).v ) - 2*(1.5*l(x,y-1,z).v-2*l(x,y-1,z-1).v+0.5*l(x,y-1,z-2).v) +1.5*(1.5*l(x,y,z).v-2*l(x,y,z-1).v+0.5*l(x,y,z-2).v);
           dwdydz1 = 0.5*( 1.5*l(x,y-2,z).w-2*l(x,y-2,z-1).w+0.5*l(x,y-2,z-2).w ) - 2*(1.5*l(x,y-1,z).w-2*l(x,y-1,z-1).w+0.5*l(x,y-1,z-2).w) +1.5*(1.5*l(x,y,z).w-2*l(x,y,z-1).w+0.5*l(x,y,z-2).w);
        }
      }

      if (normX > 0 && x+2 <    big)
      {
        dudx1 = -1.5*l(x,y,z).u+2.0*l(x+1,y,z).u-0.5*l(x+2,y,z).u;
        dvdx1 = -1.5*l(x,y,z).v+2.0*l(x+1,y,z).v-0.5*l(x+2,y,z).v;
        dwdx1 = -1.5*l(x,y,z).w+2.0*l(x+1,y,z).w-0.5*l(x+2,y,z).w;
        dudx2 =      l(x,y,z).u-2.0*l(x+1,y,z).u+    l(x+2,y,z).u;
        dvdx2 =      l(x,y,z).v-2.0*l(x+1,y,z).v+    l(x+2,y,z).v;
        dwdx2 =      l(x,y,z).w-2.0*l(x+1,y,z).w+    l(x+2,y,z).w;
      }
      if (normX < 0 && x-2 >= small)
      {
        dudx1 =  1.5*l(x,y,z).u-2.0*l(x-1,y,z).u+0.5*l(x-2,y,z).u;
        dvdx1 =  1.5*l(x,y,z).v-2.0*l(x-1,y,z).v+0.5*l(x-2,y,z).v;
        dwdx1 =  1.5*l(x,y,z).w-2.0*l(x-1,y,z).w+0.5*l(x-2,y,z).w;
        dudx2 =      l(x,y,z).u-2.0*l(x-1,y,z).u+    l(x-2,y,z).u;
        dvdx2 =      l(x,y,z).v-2.0*l(x-1,y,z).v+    l(x-2,y,z).v;
        dwdx2 =      l(x,y,z).w-2.0*l(x-1,y,z).w+    l(x-2,y,z).w;
      }
      if (normY > 0 && y+2 <    big)
      {
        dudy1 = -1.5*l(x,y,z).u+2.0*l(x,y+1,z).u-0.5*l(x,y+2,z).u;
        dvdy1 = -1.5*l(x,y,z).v+2.0*l(x,y+1,z).v-0.5*l(x,y+2,z).v;
        dwdy1 = -1.5*l(x,y,z).w+2.0*l(x,y+1,z).w-0.5*l(x,y+2,z).w;
        dudy2 =      l(x,y,z).u-2.0*l(x,y+1,z).u+    l(x,y+2,z).u;
        dvdy2 =      l(x,y,z).v-2.0*l(x,y+1,z).v+    l(x,y+2,z).v;
        dwdy2 =      l(x,y,z).w-2.0*l(x,y+1,z).w+    l(x,y+2,z).w;
      }
      if (normY < 0 && y-2 >= small)
      {
        dudy1 =  1.5*l(x,y,z).u-2.0*l(x,y-1,z).u+0.5*l(x,y-2,z).u;
        dvdy1 =  1.5*l(x,y,z).v-2.0*l(x,y-1,z).v+0.5*l(x,y-2,z).v;
        dwdy1 =  1.5*l(x,y,z).w-2.0*l(x,y-1,z).w+0.5*l(x,y-2,z).w;
        dudy2 =      l(x,y,z).u-2.0*l(x,y-1,z).u+    l(x,y-2,z).u;
        dvdy2 =      l(x,y,z).v-2.0*l(x,y-1,z).v+    l(x,y-2,z).v;
        dwdy2 =      l(x,y,z).w-2.0*l(x,y-1,z).w+    l(x,y-2,z).w;
      }
      if (normZ > 0 && z+2 <    big)
      {
        dudz1 = -1.5*l(x,y,z).u+2.0*l(x,y,z+1).u-0.5*l(x,y,z+2).u;
        dvdz1 = -1.5*l(x,y,z).v+2.0*l(x,y,z+1).v-0.5*l(x,y,z+2).v;
        dwdz1 = -1.5*l(x,y,z).w+2.0*l(x,y,z+1).w-0.5*l(x,y,z+2).w;
        dudz2 =      l(x,y,z).u-2.0*l(x,y,z+1).u+    l(x,y,z+2).u;
        dvdz2 =      l(x,y,z).v-2.0*l(x,y,z+1).v+    l(x,y,z+2).v;
        dwdz2 =      l(x,y,z).w-2.0*l(x,y,z+1).w+    l(x,y,z+2).w;
      }
      if (normZ < 0 && z-2 >= small)
      {
        dudz1 =  1.5*l(x,y,z).u-2.0*l(x,y,z-1).u+0.5*l(x,y,z-2).u;
        dvdz1 =  1.5*l(x,y,z).v-2.0*l(x,y,z-1).v+0.5*l(x,y,z-2).v;
        dwdz1 =  1.5*l(x,y,z).w-2.0*l(x,y,z-1).w+0.5*l(x,y,z-2).w;
        dudz2 =      l(x,y,z).u-2.0*l(x,y,z-1).u+    l(x,y,z-2).u;
        dvdz2 =      l(x,y,z).v-2.0*l(x,y,z-1).v+    l(x,y,z-2).v;
        dwdz2 =      l(x,y,z).w-2.0*l(x,y,z-1).w+    l(x,y,z-2).w;
      }
      if (normX > 0 && x+3 <    big)
      {
        dudx3 = -l(x,y,z).u + 3*l(x+1,y,z).u - 3*l(x+2,y,z).u + l(x+3,y,z).u; 
        dvdx3 = -l(x,y,z).v + 3*l(x+1,y,z).v - 3*l(x+2,y,z).v + l(x+3,y,z).v;
        dwdx3 = -l(x,y,z).w + 3*l(x+1,y,z).w - 3*l(x+2,y,z).w + l(x+3,y,z).w;
      }
      if (normX < 0 && x-3 >= small)
      {
        dudx3 =  l(x,y,z).u - 3*l(x-1,y,z).u + 3*l(x-2,y,z).u - l(x-3,y,z).u; 
        dvdx3 =  l(x,y,z).v - 3*l(x-1,y,z).v + 3*l(x-2,y,z).v - l(x-3,y,z).v;
        dwdx3 =  l(x,y,z).w - 3*l(x-1,y,z).w + 3*l(x-2,y,z).w - l(x-3,y,z).w;
      }
      if (normY > 0 && y+3 <    big) 
      {
        dudy3 = -l(x,y,z).u + 3*l(x,y+1,z).u - 3*l(x,y+2,z).u + l(x,y+3,z).u;
        dvdy3 = -l(x,y,z).v + 3*l(x,y+1,z).v - 3*l(x,y+2,z).v + l(x,y+3,z).v;
        dwdy3 = -l(x,y,z).w + 3*l(x,y+1,z).w - 3*l(x,y+2,z).w + l(x,y+3,z).w;
      }
      if (normY < 0 && y-3 >= small)
      {
        dudy3 =  l(x,y,z).u - 3*l(x,y-1,z).u + 3*l(x,y-2,z).u - l(x,y-3,z).u;
        dvdy3 =  l(x,y,z).v - 3*l(x,y-1,z).v + 3*l(x,y-2,z).v - l(x,y-3,z).v;
        dwdy3 =  l(x,y,z).w - 3*l(x,y-1,z).w + 3*l(x,y-2,z).w - l(x,y-3,z).w;
      }
      if (normZ > 0 && z+3 <    big) 
      {
        dudz3 = -l(x,y,z).u + 3*l(x,y,z+1).u - 3*l(x,y,z+2).u + l(x,y,z+3).u;
        dvdz3 = -l(x,y,z).v + 3*l(x,y,z+1).v - 3*l(x,y,z+2).v + l(x,y,z+3).v;
        dwdz3 = -l(x,y,z).w + 3*l(x,y,z+1).w - 3*l(x,y,z+2).w + l(x,y,z+3).w;
      }
      if (normZ < 0 && z-3 >= small)
      {
        dudz3 =  l(x,y,z).u - 3*l(x,y,z-1).u + 3*l(x,y,z-2).u - l(x,y,z-3).u;
        dvdz3 =  l(x,y,z).v - 3*l(x,y,z-1).v + 3*l(x,y,z-2).v - l(x,y,z-3).v;
        dwdz3 =  l(x,y,z).w - 3*l(x,y,z-1).w + 3*l(x,y,z-2).w - l(x,y,z-3).w;
      }

      const Real dudx = dudx1 + dudx2*(ix-x) + dudxdy1*(iy-y) + dudxdz1*(iz-z) + 0.5*dudx3*(ix-x)*(ix-x);
      const Real dvdx = dvdx1 + dvdx2*(ix-x) + dvdxdy1*(iy-y) + dvdxdz1*(iz-z) + 0.5*dvdx3*(ix-x)*(ix-x);
      const Real dwdx = dwdx1 + dwdx2*(ix-x) + dwdxdy1*(iy-y) + dwdxdz1*(iz-z) + 0.5*dwdx3*(ix-x)*(ix-x);
      const Real dudy = dudy1 + dudy2*(iy-y) + dudydz1*(iz-z) + dudxdy1*(ix-x) + 0.5*dudy3*(iy-y)*(iy-y);
      const Real dvdy = dvdy1 + dvdy2*(iy-y) + dvdydz1*(iz-z) + dvdxdy1*(ix-x) + 0.5*dvdy3*(iy-y)*(iy-y);
      const Real dwdy = dwdy1 + dwdy2*(iy-y) + dwdydz1*(iz-z) + dwdxdy1*(ix-x) + 0.5*dwdy3*(iy-y)*(iy-y);
      const Real dudz = dudz1 + dudz2*(iz-z) + dudxdz1*(ix-x) + dudydz1*(iy-y) + 0.5*dudz3*(iz-z)*(iz-z);
      const Real dvdz = dvdz1 + dvdz2*(iz-z) + dvdxdz1*(ix-x) + dvdydz1*(iy-y) + 0.5*dvdz3*(iz-z)*(iz-z);
      const Real dwdz = dwdz1 + dwdz2*(iz-z) + dwdxdz1*(ix-x) + dwdydz1*(iy-y) + 0.5*dwdz3*(iz-z)*(iz-z);       

      //normals computed with Towers 2009
      // Actually using the volume integral, since (\iint -P \hat{n} dS) = (\iiint -\nabla P dV). Also, P*\nabla\Chi = \nabla P
      // penalty-accel and surf-force match up if resolution is high enough (200 points per fish)
      const Real P = l(ix,iy,iz).p;
      //const Real fXV = D11 * normX + D12 * normY + D13 * normZ;
      //const Real fYV = D12 * normX + D22 * normY + D23 * normZ;
      //const Real fZV = D13 * normX + D23 * normY + D33 * normZ;
      const Real fXV = _1oH * (dudx * normX + dudy * normY + dudz * normZ);
      const Real fYV = _1oH * (dvdx * normX + dvdy * normY + dvdz * normZ);
      const Real fZV = _1oH * (dwdx * normX + dwdy * normY + dwdz * normZ);

      const Real fXP = -P * normX, fYP = -P * normY, fZP = -P * normZ;
      const Real fXT = fXV+fXP, fYT = fYV+fYP, fZT = fZV+fZP;

      //store:
      o->ss[i] = o->sectionMarker[iz][iy][ix];
      o->pX[i] = p[0]; o->pY[i] = p[1]; o->pZ[i] = p[2];
      o->P[i] = P; o->fX[i] = fXT; o->fY[i] = fYT; o->fZ[i] = fZT;
      o->vxDef[i] = o->udef[iz][iy][ix][0]; o->vX[i] = l(ix,iy,iz).u;
      o->vyDef[i] = o->udef[iz][iy][ix][1]; o->vY[i] = l(ix,iy,iz).v;
      o->vzDef[i] = o->udef[iz][iy][ix][2]; o->vZ[i] = l(ix,iy,iz).w;

      //additive quantities:
      // o->surface += o->surface[i]->delta;
      o->gammax += normY*o->vZ[i] - normZ*o->vY[i];
      o->gammay += normZ*o->vX[i] - normX*o->vZ[i];
      o->gammaz += normX*o->vY[i] - normY*o->vX[i];
      //forces (total, visc, pressure):
      o->forcex   += fXT; o->forcey   += fYT; o->forcez   += fZT;
      o->forcex_V += fXV; o->forcey_V += fYV; o->forcez_V += fZV;
      o->forcex_P += fXP; o->forcey_P += fYP; o->forcez_P += fZP;
      //torque:
      o->torquex  += (p[1]-CM[1])*fZT - (p[2]-CM[2])*fYT;
      o->torquey  += (p[2]-CM[2])*fXT - (p[0]-CM[0])*fZT;
      o->torquez  += (p[0]-CM[0])*fYT - (p[1]-CM[1])*fXT;
      /*
      if(tempIt->second->sectionMarker[iz][iy][ix] > 0){
        const Real * const pHinge2 = tempIt->second->hinge2LabFrame;
        (*measures)[19] += (p[1]-pHinge2[1])*fZT - (p[2]-pHinge2[2])*fYT;
        (*measures)[20] += (p[2]-pHinge2[2])*fXT - (p[0]-pHinge2[0])*fZT;
        (*measures)[21] += (p[0]-pHinge2[0])*fYT - (p[1]-pHinge2[1])*fXT;
      }
      */
      //thrust, drag:
      const Real forcePar= fXT*velUnit[0] +fYT*velUnit[1] +fZT*velUnit[2];
      o->thrust += .5*(forcePar + std::fabs(forcePar));
      o->drag   -= .5*(forcePar - std::fabs(forcePar));

      //power output (and negative definite variant which ensures no elastic energy absorption)
      // This is total power, for overcoming not only deformation, but also the oncoming velocity. Work done by fluid, not by the object (for that, just take -ve)
      const Real powOut = fXT * o->vX[i] + fYT * o->vY[i] + fZT * o->vZ[i];
      //deformation power output (and negative definite variant which ensures no elastic energy absorption)
      const Real powDef = fXT*o->vxDef[i] + fYT*o->vyDef[i] + fZT*o->vzDef[i];
      o->Pout        += powOut; o->PoutBnd     += std::min((Real)0, powOut);
      o->defPower    += powDef; o->defPowerBnd += std::min((Real)0, powDef);

      // Compute P_locomotion = Force*(uTrans + uRot)
      const Real rVec[3] = {p[0]-CM[0], p[1]-CM[1], p[2]-CM[2]};
      const Real uSolid[3] = {
	        uTrans[0] + omega[1]*rVec[2] - rVec[1]*omega[2],
	        uTrans[1] + omega[2]*rVec[0] - rVec[2]*omega[0],
	        uTrans[2] + omega[0]*rVec[1] - rVec[0]*omega[1]
      };
      o->pLocom += fXT*uSolid[0] + fYT*uSolid[1] + fZT*uSolid[2];
    }
  }
};

}

void ComputeForces::operator()(const Real dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  KernelComputeForces K(sim.nu,sim.dt,sim.obstacle_vector);
  compute<KernelComputeForces,FluidGridMPI,LabMPI,FluidGridMPI>(K,sim.grid,sim.grid);

  // do the final reductions and so on
  sim.obstacle_vector->computeForces();
}

CubismUP_3D_NAMESPACE_END

