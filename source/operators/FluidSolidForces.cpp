//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "FluidSolidForces.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

struct KernelComputeForces
{
  ObstacleVector * const obstacle_vector;
  const Real nu, dt;

  const int big   = 5;
  const int small = -4;

  const int stencil_start[3] = {small,small,small}, stencil_end[3] = {big,big,big};
  StencilInfo stencil{small,small,small, big,big,big, true, {0,1,2}};
  StencilInfo stencil2{small,small,small, big,big,big, true, {0}};
  SimulationData & sim;

  const std::vector<cubism::BlockInfo>& presInfo = sim.presInfo();

  KernelComputeForces(Real _nu, Real _dt, ObstacleVector* ov, SimulationData& s) :
    obstacle_vector(ov), nu(_nu), dt(_dt), sim(s) { }

  void operator()(VectorLab& lab, ScalarLab& chiLab, const BlockInfo& info, const BlockInfo& info2) const
  {
    for (const auto &obstacle : obstacle_vector->getObstacleVector())
      visit(lab, chiLab, info, info2, obstacle.get());
  }

  void visit(VectorLab& l, ScalarLab& chiLab, const BlockInfo& info, const BlockInfo& info2, Obstacle* const op) const
  {
    const ScalarBlock & presBlock = *(ScalarBlock*)presInfo[info.blockID].ptrBlock;
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

    //loop over elements of block info that have nonzero gradChi
    for(int i=0; i<o->nPoints; i++)
    {
      Real p[3];
      const int ix = o->surface[i]->ix;
      const int iy = o->surface[i]->iy;
      const int iz = o->surface[i]->iz;

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
        if (chiLab(x,y,z).s <0.3 && found == 0) break;

	if (ix + kk*dx + 1 >= ScalarBlock::sizeX + big-1 || ix + kk*dx -1 < small) break;
	if (iy + kk*dy + 1 >= ScalarBlock::sizeY + big-1 || iy + kk*dy -1 < small) break;
	if (iz + kk*dz + 1 >= ScalarBlock::sizeZ + big-1 || iz + kk*dz -1 < small) break;

        x  = ix + kk*dx; 
        y  = iy + kk*dy;
        z  = iz + kk*dz;
        if (chiLab(x,y,z).s < 0.3 ) {found ++; break;}
      }

      Real dudx1 = normX > 0 ? (l(x+1,y,z).u[0]-l(x,y,z).u[0]) : (l(x,y,z).u[0]-l(x-1,y,z).u[0]);
      Real dvdx1 = normX > 0 ? (l(x+1,y,z).u[1]-l(x,y,z).u[1]) : (l(x,y,z).u[1]-l(x-1,y,z).u[1]);
      Real dwdx1 = normX > 0 ? (l(x+1,y,z).u[2]-l(x,y,z).u[2]) : (l(x,y,z).u[2]-l(x-1,y,z).u[2]);
      Real dudy1 = normY > 0 ? (l(x,y+1,z).u[0]-l(x,y,z).u[0]) : (l(x,y,z).u[0]-l(x,y-1,z).u[0]);
      Real dvdy1 = normY > 0 ? (l(x,y+1,z).u[1]-l(x,y,z).u[1]) : (l(x,y,z).u[1]-l(x,y-1,z).u[1]);
      Real dwdy1 = normY > 0 ? (l(x,y+1,z).u[2]-l(x,y,z).u[2]) : (l(x,y,z).u[2]-l(x,y-1,z).u[2]);
      Real dudz1 = normY > 0 ? (l(x,y,z+1).u[0]-l(x,y,z).u[0]) : (l(x,y,z).u[0]-l(x,y,z-1).u[0]);
      Real dvdz1 = normY > 0 ? (l(x,y,z+1).u[1]-l(x,y,z).u[1]) : (l(x,y,z).u[1]-l(x,y,z-1).u[1]);
      Real dwdz1 = normY > 0 ? (l(x,y,z+1).u[2]-l(x,y,z).u[2]) : (l(x,y,z).u[2]-l(x,y,z-1).u[2]);
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

      dudxdy1 = 0.25*(l(x+1,y+1,z).u[0]+l(x-1,y-1,z).u[0]-l(x+1,y-1,z).u[0]-l(x-1,y+1,z).u[0]);
      dvdxdy1 = 0.25*(l(x+1,y+1,z).u[1]+l(x-1,y-1,z).u[1]-l(x+1,y-1,z).u[1]-l(x-1,y+1,z).u[1]);
      dwdxdy1 = 0.25*(l(x+1,y+1,z).u[2]+l(x-1,y-1,z).u[2]-l(x+1,y-1,z).u[2]-l(x-1,y+1,z).u[2]);
      if (normX > 0 && normY > 0)
      {
        dudxdy1 = (l(x+1,y+1,z).u[0]+l(x,y,z).u[0]-l(x+1,y,z).u[0]-l(x,y+1,z).u[0]);
        dvdxdy1 = (l(x+1,y+1,z).u[1]+l(x,y,z).u[1]-l(x+1,y,z).u[1]-l(x,y+1,z).u[1]);
        dwdxdy1 = (l(x+1,y+1,z).u[2]+l(x,y,z).u[2]-l(x+1,y,z).u[2]-l(x,y+1,z).u[2]);
        if ((x+2 < big) && (y+2 < big))
        {
           dudxdy1 = -0.5*( -1.5*l(x+2,y,z).u[0]+2*l(x+2,y+1,z).u[0]-0.5*l(x+2,y+2,z).u[0] ) + 2*(-1.5*l(x+1,y,z).u[0]+2*l(x+1,y+1,z).u[0]-0.5*l(x+1,y+2,z).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y+1,z).u[0]-0.5*l(x,y+2,z).u[0]);
           dvdxdy1 = -0.5*( -1.5*l(x+2,y,z).u[1]+2*l(x+2,y+1,z).u[1]-0.5*l(x+2,y+2,z).u[1] ) + 2*(-1.5*l(x+1,y,z).u[1]+2*l(x+1,y+1,z).u[1]-0.5*l(x+1,y+2,z).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y+1,z).u[1]-0.5*l(x,y+2,z).u[1]);
           dwdxdy1 = -0.5*( -1.5*l(x+2,y,z).u[2]+2*l(x+2,y+1,z).u[2]-0.5*l(x+2,y+2,z).u[2] ) + 2*(-1.5*l(x+1,y,z).u[2]+2*l(x+1,y+1,z).u[2]-0.5*l(x+1,y+2,z).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y+1,z).u[2]-0.5*l(x,y+2,z).u[2]);
        }
      }
      if (normX < 0 && normY > 0)
      {
        dudxdy1 = (l(x,y+1,z).u[0]+l(x-1,y,z).u[0]-l(x,y,z).u[0]-l(x-1,y+1,z).u[0]);
        dvdxdy1 = (l(x,y+1,z).u[1]+l(x-1,y,z).u[1]-l(x,y,z).u[1]-l(x-1,y+1,z).u[1]);
        dwdxdy1 = (l(x,y+1,z).u[2]+l(x-1,y,z).u[2]-l(x,y,z).u[2]-l(x-1,y+1,z).u[2]);
        if ((y+2 < big) && (x-2 >= small))
        {
           dudxdy1 = 0.5*( -1.5*l(x-2,y,z).u[0]+2*l(x-2,y+1,z).u[0]-0.5*l(x-2,y+2,z).u[0] ) - 2*(-1.5*l(x-1,y,z).u[0]+2*l(x-1,y+1,z).u[0]-0.5*l(x-1,y+2,z).u[0])+1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y+1,z).u[0]-0.5*l(x,y+2,z).u[0]);
           dvdxdy1 = 0.5*( -1.5*l(x-2,y,z).u[1]+2*l(x-2,y+1,z).u[1]-0.5*l(x-2,y+2,z).u[1] ) - 2*(-1.5*l(x-1,y,z).u[1]+2*l(x-1,y+1,z).u[1]-0.5*l(x-1,y+2,z).u[1])+1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y+1,z).u[1]-0.5*l(x,y+2,z).u[1]);
           dwdxdy1 = 0.5*( -1.5*l(x-2,y,z).u[2]+2*l(x-2,y+1,z).u[2]-0.5*l(x-2,y+2,z).u[2] ) - 2*(-1.5*l(x-1,y,z).u[2]+2*l(x-1,y+1,z).u[2]-0.5*l(x-1,y+2,z).u[2])+1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y+1,z).u[2]-0.5*l(x,y+2,z).u[2]);
        }
      }
      if (normX > 0 && normY < 0)
      {
        dudxdy1 = (l(x+1,y,z).u[0]+l(x,y-1,z).u[0]-l(x+1,y-1,z).u[0]-l(x,y,z).u[0]);
        dvdxdy1 = (l(x+1,y,z).u[1]+l(x,y-1,z).u[1]-l(x+1,y-1,z).u[1]-l(x,y,z).u[1]);
        dwdxdy1 = (l(x+1,y,z).u[2]+l(x,y-1,z).u[2]-l(x+1,y-1,z).u[2]-l(x,y,z).u[2]);
        if ((x+2 < big) && (y-2 >= small))
        {
           dudxdy1 = -0.5*( 1.5*l(x+2,y,z).u[0]-2*l(x+2,y-1,z).u[0]+0.5*l(x+2,y-2,z).u[0] ) + 2*(1.5*l(x+1,y,z).u[0]-2*l(x+1,y-1,z).u[0]+0.5*l(x+1,y-2,z).u[0]) -1.5*(1.5*l(x,y,z).u[0]-2*l(x,y-1,z).u[0]+0.5*l(x,y-2,z).u[0]);
           dvdxdy1 = -0.5*( 1.5*l(x+2,y,z).u[1]-2*l(x+2,y-1,z).u[1]+0.5*l(x+2,y-2,z).u[1] ) + 2*(1.5*l(x+1,y,z).u[1]-2*l(x+1,y-1,z).u[1]+0.5*l(x+1,y-2,z).u[1]) -1.5*(1.5*l(x,y,z).u[1]-2*l(x,y-1,z).u[1]+0.5*l(x,y-2,z).u[1]);
           dwdxdy1 = -0.5*( 1.5*l(x+2,y,z).u[2]-2*l(x+2,y-1,z).u[2]+0.5*l(x+2,y-2,z).u[2] ) + 2*(1.5*l(x+1,y,z).u[2]-2*l(x+1,y-1,z).u[2]+0.5*l(x+1,y-2,z).u[2]) -1.5*(1.5*l(x,y,z).u[2]-2*l(x,y-1,z).u[2]+0.5*l(x,y-2,z).u[2]);
        }
      }
      if (normX < 0 && normY < 0)
      {
        dudxdy1 = (l(x,y,z).u[0]+l(x-1,y-1,z).u[0]-l(x,y-1,z).u[0]-l(x-1,y,z).u[0]);
        dvdxdy1 = (l(x,y,z).u[1]+l(x-1,y-1,z).u[1]-l(x,y-1,z).u[1]-l(x-1,y,z).u[1]);
        dwdxdy1 = (l(x,y,z).u[2]+l(x-1,y-1,z).u[2]-l(x,y-1,z).u[2]-l(x-1,y,z).u[2]);
        if ((x-2 >= small) && (y-2 >= small))
        {
           dudxdy1 = 0.5*( 1.5*l(x-2,y,z).u[0]-2*l(x-2,y-1,z).u[0]+0.5*l(x-2,y-2,z).u[0] ) - 2*(1.5*l(x-1,y,z).u[0]-2*l(x-1,y-1,z).u[0]+0.5*l(x-1,y-2,z).u[0]) +1.5*(1.5*l(x,y,z).u[0]-2*l(x,y-1,z).u[0]+0.5*l(x,y-2,z).u[0]);
           dvdxdy1 = 0.5*( 1.5*l(x-2,y,z).u[1]-2*l(x-2,y-1,z).u[1]+0.5*l(x-2,y-2,z).u[1] ) - 2*(1.5*l(x-1,y,z).u[1]-2*l(x-1,y-1,z).u[1]+0.5*l(x-1,y-2,z).u[1]) +1.5*(1.5*l(x,y,z).u[1]-2*l(x,y-1,z).u[1]+0.5*l(x,y-2,z).u[1]);
           dwdxdy1 = 0.5*( 1.5*l(x-2,y,z).u[2]-2*l(x-2,y-1,z).u[2]+0.5*l(x-2,y-2,z).u[2] ) - 2*(1.5*l(x-1,y,z).u[2]-2*l(x-1,y-1,z).u[2]+0.5*l(x-1,y-2,z).u[2]) +1.5*(1.5*l(x,y,z).u[2]-2*l(x,y-1,z).u[2]+0.5*l(x,y-2,z).u[2]);
        }
      }

      dudxdz1 = 0.25*(l(x+1,y,z+1).u[0]+l(x-1,y,z-1).u[0]-l(x+1,y,z-1).u[0]-l(x-1,y,z+1).u[0]);
      dvdxdz1 = 0.25*(l(x+1,y,z+1).u[1]+l(x-1,y,z-1).u[1]-l(x+1,y,z-1).u[1]-l(x-1,y,z+1).u[1]);
      dwdxdz1 = 0.25*(l(x+1,y,z+1).u[2]+l(x-1,y,z-1).u[2]-l(x+1,y,z-1).u[2]-l(x-1,y,z+1).u[2]);
      if (normX > 0 && normZ > 0)
      {
        dudxdz1 = (l(x+1,y,z+1).u[0]+l(x,y,z).u[0]-l(x+1,y,z).u[0]-l(x,y,z+1).u[0]);
        dvdxdz1 = (l(x+1,y,z+1).u[1]+l(x,y,z).u[1]-l(x+1,y,z).u[1]-l(x,y,z+1).u[1]);
        dwdxdz1 = (l(x+1,y,z+1).u[2]+l(x,y,z).u[2]-l(x+1,y,z).u[2]-l(x,y,z+1).u[2]);
        if ((x+2 < big) && (z+2 < big))
        {
           dudxdz1 = -0.5*( -1.5*l(x+2,y,z).u[0]+2*l(x+2,y,z+1).u[0]-0.5*l(x+2,y,z+2).u[0] ) + 2*(-1.5*l(x+1,y,z).u[0]+2*l(x+1,y,z+1).u[0]-0.5*l(x+1,y,z+2).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+1).u[0]-0.5*l(x,y,z+2).u[0]);
           dvdxdz1 = -0.5*( -1.5*l(x+2,y,z).u[1]+2*l(x+2,y,z+1).u[1]-0.5*l(x+2,y,z+2).u[1] ) + 2*(-1.5*l(x+1,y,z).u[1]+2*l(x+1,y,z+1).u[1]-0.5*l(x+1,y,z+2).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+1).u[1]-0.5*l(x,y,z+2).u[1]);
           dwdxdz1 = -0.5*( -1.5*l(x+2,y,z).u[2]+2*l(x+2,y,z+1).u[2]-0.5*l(x+2,y,z+2).u[2] ) + 2*(-1.5*l(x+1,y,z).u[2]+2*l(x+1,y,z+1).u[2]-0.5*l(x+1,y,z+2).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+1).u[2]-0.5*l(x,y,z+2).u[2]);
        }
      }
      if (normX < 0 && normZ > 0)
      {
        dudxdz1 = (l(x,y,z+1).u[0]+l(x-1,y,z).u[0]-l(x,y,z).u[0]-l(x-1,y,z+1).u[0]);
        dvdxdz1 = (l(x,y,z+1).u[1]+l(x-1,y,z).u[1]-l(x,y,z).u[1]-l(x-1,y,z+1).u[1]);
        dwdxdz1 = (l(x,y,z+1).u[2]+l(x-1,y,z).u[2]-l(x,y,z).u[2]-l(x-1,y,z+1).u[2]);
        if ((z+2 < big) && (x-2 >= small))
        {
           dudxdz1 = 0.5*( -1.5*l(x-2,y,z).u[0]+2*l(x-2,y,z+1).u[0]-0.5*l(x-2,y,z+2).u[0] ) - 2*(-1.5*l(x-1,y,z).u[0]+2*l(x-1,y,z+1).u[0]-0.5*l(x-1,y,z+2).u[0])+1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+1).u[0]-0.5*l(x,y,z+2).u[0]);
           dvdxdz1 = 0.5*( -1.5*l(x-2,y,z).u[1]+2*l(x-2,y,z+1).u[1]-0.5*l(x-2,y,z+2).u[1] ) - 2*(-1.5*l(x-1,y,z).u[1]+2*l(x-1,y,z+1).u[1]-0.5*l(x-1,y,z+2).u[1])+1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+1).u[1]-0.5*l(x,y,z+2).u[1]);
           dwdxdz1 = 0.5*( -1.5*l(x-2,y,z).u[2]+2*l(x-2,y,z+1).u[2]-0.5*l(x-2,y,z+2).u[2] ) - 2*(-1.5*l(x-1,y,z).u[2]+2*l(x-1,y,z+1).u[2]-0.5*l(x-1,y,z+2).u[2])+1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+1).u[2]-0.5*l(x,y,z+2).u[2]);
        }
      }
      if (normX > 0 && normZ < 0)
      {
        dudxdz1 = (l(x+1,y,z).u[0]+l(x,y,z-1).u[0]-l(x+1,y,z-1).u[0]-l(x,y,z).u[0]);
        dvdxdz1 = (l(x+1,y,z).u[1]+l(x,y,z-1).u[1]-l(x+1,y,z-1).u[1]-l(x,y,z).u[1]);
        dwdxdz1 = (l(x+1,y,z).u[2]+l(x,y,z-1).u[2]-l(x+1,y,z-1).u[2]-l(x,y,z).u[2]);
        if ((x+2 < big) && (z-2 >= small))
        {
           dudxdz1 = -0.5*( 1.5*l(x+2,y,z).u[0]-2*l(x+2,y,z-1).u[0]+0.5*l(x+2,y,z-2).u[0] ) + 2*(1.5*l(x+1,y,z).u[0]-2*l(x+1,y,z-1).u[0]+0.5*l(x+1,y,z-2).u[0]) -1.5*(1.5*l(x,y,z).u[0]-2*l(x,y,z-1).u[0]+0.5*l(x,y,z-2).u[0]);
           dvdxdz1 = -0.5*( 1.5*l(x+2,y,z).u[1]-2*l(x+2,y,z-1).u[1]+0.5*l(x+2,y,z-2).u[1] ) + 2*(1.5*l(x+1,y,z).u[1]-2*l(x+1,y,z-1).u[1]+0.5*l(x+1,y,z-2).u[1]) -1.5*(1.5*l(x,y,z).u[1]-2*l(x,y,z-1).u[1]+0.5*l(x,y,z-2).u[1]);
           dwdxdz1 = -0.5*( 1.5*l(x+2,y,z).u[2]-2*l(x+2,y,z-1).u[2]+0.5*l(x+2,y,z-2).u[2] ) + 2*(1.5*l(x+1,y,z).u[2]-2*l(x+1,y,z-1).u[2]+0.5*l(x+1,y,z-2).u[2]) -1.5*(1.5*l(x,y,z).u[2]-2*l(x,y,z-1).u[2]+0.5*l(x,y,z-2).u[2]);
        }
      }
      if (normX < 0 && normZ < 0)
      {
        dudxdz1 = (l(x,y,z).u[0]+l(x-1,y,z-1).u[0]-l(x,y,z-1).u[0]-l(x-1,y,z).u[0]);
        dvdxdz1 = (l(x,y,z).u[1]+l(x-1,y,z-1).u[1]-l(x,y,z-1).u[1]-l(x-1,y,z).u[1]);
        dwdxdz1 = (l(x,y,z).u[2]+l(x-1,y,z-1).u[2]-l(x,y,z-1).u[2]-l(x-1,y,z).u[2]);
        if ((x-2 >= small) && (z-2 >= small))
        {
           dudxdz1 = 0.5*( 1.5*l(x-2,y,z).u[0]-2*l(x-2,y,z-1).u[0]+0.5*l(x-2,y,z-2).u[0] ) - 2*(1.5*l(x-1,y,z).u[0]-2*l(x-1,y,z-1).u[0]+0.5*l(x-1,y,z-2).u[0]) +1.5*(1.5*l(x,y,z).u[0]-2*l(x,y,z-1).u[0]+0.5*l(x,y,z-2).u[0]);
           dvdxdz1 = 0.5*( 1.5*l(x-2,y,z).u[1]-2*l(x-2,y,z-1).u[1]+0.5*l(x-2,y,z-2).u[1] ) - 2*(1.5*l(x-1,y,z).u[1]-2*l(x-1,y,z-1).u[1]+0.5*l(x-1,y,z-2).u[1]) +1.5*(1.5*l(x,y,z).u[1]-2*l(x,y,z-1).u[1]+0.5*l(x,y,z-2).u[1]);
           dwdxdz1 = 0.5*( 1.5*l(x-2,y,z).u[2]-2*l(x-2,y,z-1).u[2]+0.5*l(x-2,y,z-2).u[2] ) - 2*(1.5*l(x-1,y,z).u[2]-2*l(x-1,y,z-1).u[2]+0.5*l(x-1,y,z-2).u[2]) +1.5*(1.5*l(x,y,z).u[2]-2*l(x,y,z-1).u[2]+0.5*l(x,y,z-2).u[2]);
        }
      }

      dudydz1 = 0.25*(l(x,y+1,z+1).u[0]+l(x,y-1,z-1).u[0]-l(x,y+1,z-1).u[0]-l(x,y-1,z+1).u[0]);
      dvdydz1 = 0.25*(l(x,y+1,z+1).u[1]+l(x,y-1,z-1).u[1]-l(x,y+1,z-1).u[1]-l(x,y-1,z+1).u[1]);
      dwdydz1 = 0.25*(l(x,y+1,z+1).u[2]+l(x,y-1,z-1).u[2]-l(x,y+1,z-1).u[2]-l(x,y-1,z+1).u[2]);
      if (normY > 0 && normZ > 0)
      {
        dudydz1 = (l(x,y+1,z+1).u[0]+l(x,y,z).u[0]-l(x,y+1,z).u[0]-l(x,y,z+1).u[0]);
        dvdydz1 = (l(x,y+1,z+1).u[1]+l(x,y,z).u[1]-l(x,y+1,z).u[1]-l(x,y,z+1).u[1]);
        dwdydz1 = (l(x,y+1,z+1).u[2]+l(x,y,z).u[2]-l(x,y+1,z).u[2]-l(x,y,z+1).u[2]);
        if ((y+2 < big) && (z+2 < big))
        {
           dudydz1 = -0.5*( -1.5*l(x,y+2,z).u[0]+2*l(x,y+2,z+1).u[0]-0.5*l(x,y+2,z+2).u[0] ) + 2*(-1.5*l(x,y+1,z).u[0]+2*l(x,y+1,z+1).u[0]-0.5*l(x,y+1,z+2).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+1).u[0]-0.5*l(x,y,z+2).u[0]);
           dvdydz1 = -0.5*( -1.5*l(x,y+2,z).u[1]+2*l(x,y+2,z+1).u[1]-0.5*l(x,y+2,z+2).u[1] ) + 2*(-1.5*l(x,y+1,z).u[1]+2*l(x,y+1,z+1).u[1]-0.5*l(x,y+1,z+2).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+1).u[1]-0.5*l(x,y,z+2).u[1]);
           dwdydz1 = -0.5*( -1.5*l(x,y+2,z).u[2]+2*l(x,y+2,z+1).u[2]-0.5*l(x,y+2,z+2).u[2] ) + 2*(-1.5*l(x,y+1,z).u[2]+2*l(x,y+1,z+1).u[2]-0.5*l(x,y+1,z+2).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+1).u[2]-0.5*l(x,y,z+2).u[2]);
        }
      }
      if (normY < 0 && normZ > 0)
      {
        dudydz1 = (l(x,y,z+1).u[0]+l(x,y-1,z).u[0]-l(x,y,z).u[0]-l(x,y-1,z+1).u[0]);
        dvdydz1 = (l(x,y,z+1).u[1]+l(x,y-1,z).u[1]-l(x,y,z).u[1]-l(x,y-1,z+1).u[1]);
        dwdydz1 = (l(x,y,z+1).u[2]+l(x,y-1,z).u[2]-l(x,y,z).u[2]-l(x,y-1,z+1).u[2]);
        if ((z+2 < big) && (y-2 >= small))
        {
           dudydz1 = 0.5*( -1.5*l(x,y-2,z).u[0]+2*l(x,y-2,z+1).u[0]-0.5*l(x,y-2,z+2).u[0] ) - 2*(-1.5*l(x,y-1,z).u[0]+2*l(x,y-1,z+1).u[0]-0.5*l(x,y-1,z+2).u[0])+1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+1).u[0]-0.5*l(x,y,z+2).u[0]);
           dvdydz1 = 0.5*( -1.5*l(x,y-2,z).u[1]+2*l(x,y-2,z+1).u[1]-0.5*l(x,y-2,z+2).u[1] ) - 2*(-1.5*l(x,y-1,z).u[1]+2*l(x,y-1,z+1).u[1]-0.5*l(x,y-1,z+2).u[1])+1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+1).u[1]-0.5*l(x,y,z+2).u[1]);
           dwdydz1 = 0.5*( -1.5*l(x,y-2,z).u[2]+2*l(x,y-2,z+1).u[2]-0.5*l(x,y-2,z+2).u[2] ) - 2*(-1.5*l(x,y-1,z).u[2]+2*l(x,y-1,z+1).u[2]-0.5*l(x,y-1,z+2).u[2])+1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+1).u[2]-0.5*l(x,y,z+2).u[2]);
        }
      }
      if (normY > 0 && normZ < 0)
      {
        dudydz1 = (l(x,y+1,z).u[0]+l(x,y,z-1).u[0]-l(x,y+1,z-1).u[0]-l(x,y,z).u[0]);
        dvdydz1 = (l(x,y+1,z).u[1]+l(x,y,z-1).u[1]-l(x,y+1,z-1).u[1]-l(x,y,z).u[1]);
        dwdydz1 = (l(x,y+1,z).u[2]+l(x,y,z-1).u[2]-l(x,y+1,z-1).u[2]-l(x,y,z).u[2]);
        if ((y+2 < big) && (z-2 >= small))
        {
           dudydz1 = -0.5*( 1.5*l(x,y+2,z).u[0]-2*l(x,y+2,z-1).u[0]+0.5*l(x,y+2,z-2).u[0] ) + 2*(1.5*l(x,y+1,z).u[0]-2*l(x,y+1,z-1).u[0]+0.5*l(x,y+1,z-2).u[0]) -1.5*(1.5*l(x,y,z).u[0]-2*l(x,y,z-1).u[0]+0.5*l(x,y,z-2).u[0]);
           dvdydz1 = -0.5*( 1.5*l(x,y+2,z).u[1]-2*l(x,y+2,z-1).u[1]+0.5*l(x,y+2,z-2).u[1] ) + 2*(1.5*l(x,y+1,z).u[1]-2*l(x,y+1,z-1).u[1]+0.5*l(x,y+1,z-2).u[1]) -1.5*(1.5*l(x,y,z).u[1]-2*l(x,y,z-1).u[1]+0.5*l(x,y,z-2).u[1]);
           dwdydz1 = -0.5*( 1.5*l(x,y+2,z).u[2]-2*l(x,y+2,z-1).u[2]+0.5*l(x,y+2,z-2).u[2] ) + 2*(1.5*l(x,y+1,z).u[2]-2*l(x,y+1,z-1).u[2]+0.5*l(x,y+1,z-2).u[2]) -1.5*(1.5*l(x,y,z).u[2]-2*l(x,y,z-1).u[2]+0.5*l(x,y,z-2).u[2]);
        }
      }
      if (normY < 0 && normZ < 0)
      {
        dudydz1 = (l(x,y,z).u[0]+l(x,y-1,z-1).u[0]-l(x,y,z-1).u[0]-l(x,y-1,z).u[0]);
        dvdydz1 = (l(x,y,z).u[1]+l(x,y-1,z-1).u[1]-l(x,y,z-1).u[1]-l(x,y-1,z).u[1]);
        dwdydz1 = (l(x,y,z).u[2]+l(x,y-1,z-1).u[2]-l(x,y,z-1).u[2]-l(x,y-1,z).u[2]);
        if ((y-2 >= small) && (z-2 >= small))
        {
           dudydz1 = 0.5*( 1.5*l(x,y-2,z).u[0]-2*l(x,y-2,z-1).u[0]+0.5*l(x,y-2,z-2).u[0] ) - 2*(1.5*l(x,y-1,z).u[0]-2*l(x,y-1,z-1).u[0]+0.5*l(x,y-1,z-2).u[0]) +1.5*(1.5*l(x,y,z).u[0]-2*l(x,y,z-1).u[0]+0.5*l(x,y,z-2).u[0]);
           dvdydz1 = 0.5*( 1.5*l(x,y-2,z).u[1]-2*l(x,y-2,z-1).u[1]+0.5*l(x,y-2,z-2).u[1] ) - 2*(1.5*l(x,y-1,z).u[1]-2*l(x,y-1,z-1).u[1]+0.5*l(x,y-1,z-2).u[1]) +1.5*(1.5*l(x,y,z).u[1]-2*l(x,y,z-1).u[1]+0.5*l(x,y,z-2).u[1]);
           dwdydz1 = 0.5*( 1.5*l(x,y-2,z).u[2]-2*l(x,y-2,z-1).u[2]+0.5*l(x,y-2,z-2).u[2] ) - 2*(1.5*l(x,y-1,z).u[2]-2*l(x,y-1,z-1).u[2]+0.5*l(x,y-1,z-2).u[2]) +1.5*(1.5*l(x,y,z).u[2]-2*l(x,y,z-1).u[2]+0.5*l(x,y,z-2).u[2]);
        }
      }

      if (normX > 0 && x+2 <    big)
      {
        dudx1 = -1.5*l(x,y,z).u[0]+2.0*l(x+1,y,z).u[0]-0.5*l(x+2,y,z).u[0];
        dvdx1 = -1.5*l(x,y,z).u[1]+2.0*l(x+1,y,z).u[1]-0.5*l(x+2,y,z).u[1];
        dwdx1 = -1.5*l(x,y,z).u[2]+2.0*l(x+1,y,z).u[2]-0.5*l(x+2,y,z).u[2];
        dudx2 =      l(x,y,z).u[0]-2.0*l(x+1,y,z).u[0]+    l(x+2,y,z).u[0];
        dvdx2 =      l(x,y,z).u[1]-2.0*l(x+1,y,z).u[1]+    l(x+2,y,z).u[1];
        dwdx2 =      l(x,y,z).u[2]-2.0*l(x+1,y,z).u[2]+    l(x+2,y,z).u[2];
      }
      if (normX < 0 && x-2 >= small)
      {
        dudx1 =  1.5*l(x,y,z).u[0]-2.0*l(x-1,y,z).u[0]+0.5*l(x-2,y,z).u[0];
        dvdx1 =  1.5*l(x,y,z).u[1]-2.0*l(x-1,y,z).u[1]+0.5*l(x-2,y,z).u[1];
        dwdx1 =  1.5*l(x,y,z).u[2]-2.0*l(x-1,y,z).u[2]+0.5*l(x-2,y,z).u[2];
        dudx2 =      l(x,y,z).u[0]-2.0*l(x-1,y,z).u[0]+    l(x-2,y,z).u[0];
        dvdx2 =      l(x,y,z).u[1]-2.0*l(x-1,y,z).u[1]+    l(x-2,y,z).u[1];
        dwdx2 =      l(x,y,z).u[2]-2.0*l(x-1,y,z).u[2]+    l(x-2,y,z).u[2];
      }
      if (normY > 0 && y+2 <    big)
      {
        dudy1 = -1.5*l(x,y,z).u[0]+2.0*l(x,y+1,z).u[0]-0.5*l(x,y+2,z).u[0];
        dvdy1 = -1.5*l(x,y,z).u[1]+2.0*l(x,y+1,z).u[1]-0.5*l(x,y+2,z).u[1];
        dwdy1 = -1.5*l(x,y,z).u[2]+2.0*l(x,y+1,z).u[2]-0.5*l(x,y+2,z).u[2];
        dudy2 =      l(x,y,z).u[0]-2.0*l(x,y+1,z).u[0]+    l(x,y+2,z).u[0];
        dvdy2 =      l(x,y,z).u[1]-2.0*l(x,y+1,z).u[1]+    l(x,y+2,z).u[1];
        dwdy2 =      l(x,y,z).u[2]-2.0*l(x,y+1,z).u[2]+    l(x,y+2,z).u[2];
      }
      if (normY < 0 && y-2 >= small)
      {
        dudy1 =  1.5*l(x,y,z).u[0]-2.0*l(x,y-1,z).u[0]+0.5*l(x,y-2,z).u[0];
        dvdy1 =  1.5*l(x,y,z).u[1]-2.0*l(x,y-1,z).u[1]+0.5*l(x,y-2,z).u[1];
        dwdy1 =  1.5*l(x,y,z).u[2]-2.0*l(x,y-1,z).u[2]+0.5*l(x,y-2,z).u[2];
        dudy2 =      l(x,y,z).u[0]-2.0*l(x,y-1,z).u[0]+    l(x,y-2,z).u[0];
        dvdy2 =      l(x,y,z).u[1]-2.0*l(x,y-1,z).u[1]+    l(x,y-2,z).u[1];
        dwdy2 =      l(x,y,z).u[2]-2.0*l(x,y-1,z).u[2]+    l(x,y-2,z).u[2];
      }
      if (normZ > 0 && z+2 <    big)
      {
        dudz1 = -1.5*l(x,y,z).u[0]+2.0*l(x,y,z+1).u[0]-0.5*l(x,y,z+2).u[0];
        dvdz1 = -1.5*l(x,y,z).u[1]+2.0*l(x,y,z+1).u[1]-0.5*l(x,y,z+2).u[1];
        dwdz1 = -1.5*l(x,y,z).u[2]+2.0*l(x,y,z+1).u[2]-0.5*l(x,y,z+2).u[2];
        dudz2 =      l(x,y,z).u[0]-2.0*l(x,y,z+1).u[0]+    l(x,y,z+2).u[0];
        dvdz2 =      l(x,y,z).u[1]-2.0*l(x,y,z+1).u[1]+    l(x,y,z+2).u[1];
        dwdz2 =      l(x,y,z).u[2]-2.0*l(x,y,z+1).u[2]+    l(x,y,z+2).u[2];
      }
      if (normZ < 0 && z-2 >= small)
      {
        dudz1 =  1.5*l(x,y,z).u[0]-2.0*l(x,y,z-1).u[0]+0.5*l(x,y,z-2).u[0];
        dvdz1 =  1.5*l(x,y,z).u[1]-2.0*l(x,y,z-1).u[1]+0.5*l(x,y,z-2).u[1];
        dwdz1 =  1.5*l(x,y,z).u[2]-2.0*l(x,y,z-1).u[2]+0.5*l(x,y,z-2).u[2];
        dudz2 =      l(x,y,z).u[0]-2.0*l(x,y,z-1).u[0]+    l(x,y,z-2).u[0];
        dvdz2 =      l(x,y,z).u[1]-2.0*l(x,y,z-1).u[1]+    l(x,y,z-2).u[1];
        dwdz2 =      l(x,y,z).u[2]-2.0*l(x,y,z-1).u[2]+    l(x,y,z-2).u[2];
      }
      if (normX > 0 && x+3 <    big)
      {
        dudx3 = -l(x,y,z).u[0] + 3*l(x+1,y,z).u[0] - 3*l(x+2,y,z).u[0] + l(x+3,y,z).u[0]; 
        dvdx3 = -l(x,y,z).u[1] + 3*l(x+1,y,z).u[1] - 3*l(x+2,y,z).u[1] + l(x+3,y,z).u[1];
        dwdx3 = -l(x,y,z).u[2] + 3*l(x+1,y,z).u[2] - 3*l(x+2,y,z).u[2] + l(x+3,y,z).u[2];
      }
      if (normX < 0 && x-3 >= small)
      {
        dudx3 =  l(x,y,z).u[0] - 3*l(x-1,y,z).u[0] + 3*l(x-2,y,z).u[0] - l(x-3,y,z).u[0]; 
        dvdx3 =  l(x,y,z).u[1] - 3*l(x-1,y,z).u[1] + 3*l(x-2,y,z).u[1] - l(x-3,y,z).u[1];
        dwdx3 =  l(x,y,z).u[2] - 3*l(x-1,y,z).u[2] + 3*l(x-2,y,z).u[2] - l(x-3,y,z).u[2];
      }
      if (normY > 0 && y+3 <    big) 
      {
        dudy3 = -l(x,y,z).u[0] + 3*l(x,y+1,z).u[0] - 3*l(x,y+2,z).u[0] + l(x,y+3,z).u[0];
        dvdy3 = -l(x,y,z).u[1] + 3*l(x,y+1,z).u[1] - 3*l(x,y+2,z).u[1] + l(x,y+3,z).u[1];
        dwdy3 = -l(x,y,z).u[2] + 3*l(x,y+1,z).u[2] - 3*l(x,y+2,z).u[2] + l(x,y+3,z).u[2];
      }
      if (normY < 0 && y-3 >= small)
      {
        dudy3 =  l(x,y,z).u[0] - 3*l(x,y-1,z).u[0] + 3*l(x,y-2,z).u[0] - l(x,y-3,z).u[0];
        dvdy3 =  l(x,y,z).u[1] - 3*l(x,y-1,z).u[1] + 3*l(x,y-2,z).u[1] - l(x,y-3,z).u[1];
        dwdy3 =  l(x,y,z).u[2] - 3*l(x,y-1,z).u[2] + 3*l(x,y-2,z).u[2] - l(x,y-3,z).u[2];
      }
      if (normZ > 0 && z+3 <    big) 
      {
        dudz3 = -l(x,y,z).u[0] + 3*l(x,y,z+1).u[0] - 3*l(x,y,z+2).u[0] + l(x,y,z+3).u[0];
        dvdz3 = -l(x,y,z).u[1] + 3*l(x,y,z+1).u[1] - 3*l(x,y,z+2).u[1] + l(x,y,z+3).u[1];
        dwdz3 = -l(x,y,z).u[2] + 3*l(x,y,z+1).u[2] - 3*l(x,y,z+2).u[2] + l(x,y,z+3).u[2];
      }
      if (normZ < 0 && z-3 >= small)
      {
        dudz3 =  l(x,y,z).u[0] - 3*l(x,y,z-1).u[0] + 3*l(x,y,z-2).u[0] - l(x,y,z-3).u[0];
        dvdz3 =  l(x,y,z).u[1] - 3*l(x,y,z-1).u[1] + 3*l(x,y,z-2).u[1] - l(x,y,z-3).u[1];
        dwdz3 =  l(x,y,z).u[2] - 3*l(x,y,z-1).u[2] + 3*l(x,y,z-2).u[2] - l(x,y,z-3).u[2];
      }

      const Real dudx = dudx1 + dudx2*(ix-x) + dudxdy1*(iy-y) + dudxdz1*(iz-z);
      const Real dvdx = dvdx1 + dvdx2*(ix-x) + dvdxdy1*(iy-y) + dvdxdz1*(iz-z);
      const Real dwdx = dwdx1 + dwdx2*(ix-x) + dwdxdy1*(iy-y) + dwdxdz1*(iz-z);
      const Real dudy = dudy1 + dudy2*(iy-y) + dudydz1*(iz-z) + dudxdy1*(ix-x);
      const Real dvdy = dvdy1 + dvdy2*(iy-y) + dvdydz1*(iz-z) + dvdxdy1*(ix-x);
      const Real dwdy = dwdy1 + dwdy2*(iy-y) + dwdydz1*(iz-z) + dwdxdy1*(ix-x);
      const Real dudz = dudz1 + dudz2*(iz-z) + dudxdz1*(ix-x) + dudydz1*(iy-y);
      const Real dvdz = dvdz1 + dvdz2*(iz-z) + dvdxdz1*(ix-x) + dvdydz1*(iy-y);
      const Real dwdz = dwdz1 + dwdz2*(iz-z) + dwdxdz1*(ix-x) + dwdydz1*(iy-y);

      //normals computed with Towers 2009
      // Actually using the volume integral, since (\iint -P \hat{n} dS) = (\iiint -\nabla P dV). Also, P*\nabla\Chi = \nabla P
      // penalty-accel and surf-force match up if resolution is high enough (200 points per fish)
      const Real P = presBlock(ix,iy,iz).s;
      const Real fXV = _1oH * (dudx * normX + dudy * normY + dudz * normZ);
      const Real fYV = _1oH * (dvdx * normX + dvdy * normY + dvdz * normZ);
      const Real fZV = _1oH * (dwdx * normX + dwdy * normY + dwdz * normZ);

      const Real fXP = -P * normX, fYP = -P * normY, fZP = -P * normZ;
      const Real fXT = fXV+fXP, fYT = fYV+fYP, fZT = fZV+fZP;

      //store:
      o->pX[i] = p[0]; o->pY[i] = p[1]; o->pZ[i] = p[2];
      o->P[i] = P; o->fX[i] = fXT; o->fY[i] = fYT; o->fZ[i] = fZT;
      o->fxV[i] = _1oH * (dudx * dx + dudy * dy + dudz * dz); 
      o->fyV[i] = _1oH * (dvdx * dx + dvdy * dy + dvdz * dz); 
      o->fzV[i] = _1oH * (dwdx * dx + dwdy * dy + dwdz * dz); 

      o->vxDef[i] = o->udef[iz][iy][ix][0]; o->vX[i] = l(ix,iy,iz).u[0];
      o->vyDef[i] = o->udef[iz][iy][ix][1]; o->vY[i] = l(ix,iy,iz).u[1];
      o->vzDef[i] = o->udef[iz][iy][ix][2]; o->vZ[i] = l(ix,iy,iz).u[2];

      //forces (total, visc, pressure):
      o->forcex   += fXT; o->forcey   += fYT; o->forcez   += fZT;
      o->forcex_V += fXV; o->forcey_V += fYV; o->forcez_V += fZV;
      o->forcex_P += fXP; o->forcey_P += fYP; o->forcez_P += fZP;
      //torque:
      o->torquex  += (p[1]-CM[1])*fZT - (p[2]-CM[2])*fYT;
      o->torquey  += (p[2]-CM[2])*fXT - (p[0]-CM[0])*fZT;
      o->torquez  += (p[0]-CM[0])*fYT - (p[1]-CM[1])*fXT;
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

  KernelComputeForces K(sim.nu,sim.dt,sim.obstacle_vector,sim);
  cubism::compute<KernelComputeForces,VectorGrid,VectorLab,ScalarGrid,ScalarLab>(K,*sim.vel,*sim.chi);
  // do the final reductions and so on
  sim.obstacle_vector->computeForces();
}

CubismUP_3D_NAMESPACE_END

