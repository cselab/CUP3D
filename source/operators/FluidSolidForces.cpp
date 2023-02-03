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
  const int big   = 5;
  const int small = -4;
  const int bigg = ScalarBlock::sizeX + big-1;
  const int stencil_start[3] = {small,small,small}, stencil_end[3] = {big,big,big};
  const Real c0 = -137./60.;
  const Real c1 =    5.    ;
  const Real c2 = -  5.    ;
  const Real c3 =   10./ 3.;
  const Real c4 = -  5./ 4.;
  const Real c5 =    1./ 5.;

  inline bool inrange(const int i) const
  {
    return (i >= small && i < bigg);
  }

  StencilInfo stencil{small,small,small, big,big,big, true, {0,1,2}};
  StencilInfo stencil2{small,small,small, big,big,big, true, {0}};
  SimulationData & sim;

  const std::vector<cubism::BlockInfo>& presInfo = sim.presInfo();

  KernelComputeForces(SimulationData& s) : sim(s) {}

  void operator()(VectorLab& lab, ScalarLab& chiLab, const BlockInfo& info, const BlockInfo& info2) const
  {
    for (const auto &obstacle : sim.obstacle_vector->getObstacleVector())
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
    const Real vel_norm = std::sqrt(uTrans[0]*uTrans[0] + uTrans[1]*uTrans[1] + uTrans[2]*uTrans[2]);
    if (vel_norm>1e-9)
    {
      velUnit[0] = uTrans[0] / vel_norm;
      velUnit[1] = uTrans[1] / vel_norm;
      velUnit[2] = uTrans[2] / vel_norm;
    }

    const Real _1oH = sim.nu / info.h;

    //loop over elements of block info that have nonzero gradChi
    for(int i=0; i<o->nPoints; i++)
    {
      const int ix = o->surface[i]->ix;
      const int iy = o->surface[i]->iy;
      const int iz = o->surface[i]->iz;

      Real p[3];
      info.pos(p, ix, iy, iz);

      //shear stresses
      const Real normX = o->surface[i]->dchidx; //*h^3 (multiplied in dchidx)
      const Real normY = o->surface[i]->dchidy; //*h^3 (multiplied in dchidy)
      const Real normZ = o->surface[i]->dchidz; //*h^3 (multiplied in dchidz)
      const Real norm = 1.0/std::sqrt(normX*normX+normY*normY+normZ*normZ);
      const Real dx = normX*norm;
      const Real dy = normY*norm;
      const Real dz = normZ*norm;

      int x = ix;
      int y = iy;
      int z = iz;
      for (int kk = 0 ; kk < 5 ; kk++) //5 is arbitrary
      {
        const int dxi = round(kk*dx);
        const int dyi = round(kk*dy);
        const int dzi = round(kk*dz);
	      if (ix + dxi + 1 >= ScalarBlock::sizeX + big-1 || ix + dxi -1 < small) continue;
	      if (iy + dyi + 1 >= ScalarBlock::sizeY + big-1 || iy + dyi -1 < small) continue;
	      if (iz + dzi + 1 >= ScalarBlock::sizeZ + big-1 || iz + dzi -1 < small) continue;
        x  = ix + dxi; 
        y  = iy + dyi;
        z  = iz + dzi;
        if (chiLab(x,y,z).s < 0.01 ) break;
      }

      Real dudx1,dvdx1,dwdx1;
      Real dudy1,dvdy1,dwdy1;
      Real dudz1,dvdz1,dwdz1;      
      const int sx = normX > 0 ? +1:-1;
      const int sy = normY > 0 ? +1:-1;
      const int sz = normZ > 0 ? +1:-1;
      if (inrange(x+5*sx))
      {
        dudx1 = sx*(c0*l(x,y,z).u[0]+c1*l(x+sx,y,z).u[0]+c2*l(x+2*sx,y,z).u[0]+c3*l(x+3*sx,y,z).u[0]+c4*l(x+4*sx,y,z).u[0]+c5*l(x+5*sx,y,z).u[0]);
        dvdx1 = sx*(c0*l(x,y,z).u[1]+c1*l(x+sx,y,z).u[1]+c2*l(x+2*sx,y,z).u[1]+c3*l(x+3*sx,y,z).u[1]+c4*l(x+4*sx,y,z).u[1]+c5*l(x+5*sx,y,z).u[1]);
        dwdx1 = sx*(c0*l(x,y,z).u[2]+c1*l(x+sx,y,z).u[2]+c2*l(x+2*sx,y,z).u[2]+c3*l(x+3*sx,y,z).u[2]+c4*l(x+4*sx,y,z).u[2]+c5*l(x+5*sx,y,z).u[2]);        
      }
      else if (inrange(x+2*sx))
      {
        dudx1 = sx*(-1.5*l(x,y,z).u[0]+2.0*l(x+sx,y,z).u[0]-0.5*l(x+2*sx,y,z).u[0]);
        dvdx1 = sx*(-1.5*l(x,y,z).u[1]+2.0*l(x+sx,y,z).u[1]-0.5*l(x+2*sx,y,z).u[1]);
        dwdx1 = sx*(-1.5*l(x,y,z).u[2]+2.0*l(x+sx,y,z).u[2]-0.5*l(x+2*sx,y,z).u[2]);
      }
      else
      {
        dudx1 = sx*(l(x+sx,y,z).u[0]-l(x,y,z).u[0]);
        dvdx1 = sx*(l(x+sx,y,z).u[1]-l(x,y,z).u[1]);
        dwdx1 = sx*(l(x+sx,y,z).u[2]-l(x,y,z).u[2]);
      }
      if (inrange(y+5*sy))
      {
        dudy1 = sy*(c0*l(x,y,z).u[0]+c1*l(x,y+sy,z).u[0]+c2*l(x,y+2*sy,z).u[0]+c3*l(x,y+3*sy,z).u[0]+c4*l(x,y+4*sy,z).u[0]+c5*l(x,y+5*sy,z).u[0]);
        dvdy1 = sy*(c0*l(x,y,z).u[1]+c1*l(x,y+sy,z).u[1]+c2*l(x,y+2*sy,z).u[1]+c3*l(x,y+3*sy,z).u[1]+c4*l(x,y+4*sy,z).u[1]+c5*l(x,y+5*sy,z).u[1]);
        dwdy1 = sy*(c0*l(x,y,z).u[2]+c1*l(x,y+sy,z).u[2]+c2*l(x,y+2*sy,z).u[2]+c3*l(x,y+3*sy,z).u[2]+c4*l(x,y+4*sy,z).u[2]+c5*l(x,y+5*sy,z).u[2]);        
      }
      else if (inrange(y+2*sy))
      {
        dudy1 = sy*(-1.5*l(x,y,z).u[0]+2.0*l(x,y+sy,z).u[0]-0.5*l(x,y+2*sy,z).u[0]);
        dvdy1 = sy*(-1.5*l(x,y,z).u[1]+2.0*l(x,y+sy,z).u[1]-0.5*l(x,y+2*sy,z).u[1]);
        dwdy1 = sy*(-1.5*l(x,y,z).u[2]+2.0*l(x,y+sy,z).u[2]-0.5*l(x,y+2*sy,z).u[2]);
      }
      else
      {
        dudy1 = sx*(l(x,y+sy,z).u[0]-l(x,y,z).u[0]);
        dvdy1 = sx*(l(x,y+sy,z).u[1]-l(x,y,z).u[1]);
        dwdy1 = sx*(l(x,y+sy,z).u[2]-l(x,y,z).u[2]);
      }
      if (inrange(z+5*sz))
      {
        dudz1 = sz*(c0*l(x,y,z).u[0]+c1*l(x,y,z+sz).u[0]+c2*l(x,y,z+2*sz).u[0]+c3*l(x,y,z+3*sz).u[0]+c4*l(x,y,z+4*sz).u[0]+c5*l(x,y,z+5*sz).u[0]);
        dvdz1 = sz*(c0*l(x,y,z).u[1]+c1*l(x,y,z+sz).u[1]+c2*l(x,y,z+2*sz).u[1]+c3*l(x,y,z+3*sz).u[1]+c4*l(x,y,z+4*sz).u[1]+c5*l(x,y,z+5*sz).u[1]);
        dwdz1 = sz*(c0*l(x,y,z).u[2]+c1*l(x,y,z+sz).u[2]+c2*l(x,y,z+2*sz).u[2]+c3*l(x,y,z+3*sz).u[2]+c4*l(x,y,z+4*sz).u[2]+c5*l(x,y,z+5*sz).u[2]);        
      }
      else if (inrange(z+2*sz))
      {
        dudz1 = sz*(-1.5*l(x,y,z).u[0]+2.0*l(x,y,z+sz).u[0]-0.5*l(x,y,z+2*sz).u[0]);
        dvdz1 = sz*(-1.5*l(x,y,z).u[1]+2.0*l(x,y,z+sz).u[1]-0.5*l(x,y,z+2*sz).u[1]);
        dwdz1 = sz*(-1.5*l(x,y,z).u[2]+2.0*l(x,y,z+sz).u[2]-0.5*l(x,y,z+2*sz).u[2]);
      }
      else
      {
        dudz1 = sz*(l(x,y,z+sz).u[0]-l(x,y,z).u[0]);
        dvdz1 = sz*(l(x,y,z+sz).u[1]-l(x,y,z).u[1]);
        dwdz1 = sz*(l(x,y,z+sz).u[2]-l(x,y,z).u[2]);
      }

      const Real dudx2 = l(x-1,y,z).u[0]-2.0*l(x,y,z).u[0]+ l(x+1,y,z).u[0];
      const Real dvdx2 = l(x-1,y,z).u[1]-2.0*l(x,y,z).u[1]+ l(x+1,y,z).u[1];
      const Real dwdx2 = l(x-1,y,z).u[2]-2.0*l(x,y,z).u[2]+ l(x+1,y,z).u[2];
      const Real dudy2 = l(x,y-1,z).u[0]-2.0*l(x,y,z).u[0]+ l(x,y+1,z).u[0];
      const Real dvdy2 = l(x,y-1,z).u[1]-2.0*l(x,y,z).u[1]+ l(x,y+1,z).u[1];
      const Real dwdy2 = l(x,y-1,z).u[2]-2.0*l(x,y,z).u[2]+ l(x,y+1,z).u[2];
      const Real dudz2 = l(x,y,z-1).u[0]-2.0*l(x,y,z).u[0]+ l(x,y,z+1).u[0];
      const Real dvdz2 = l(x,y,z-1).u[1]-2.0*l(x,y,z).u[1]+ l(x,y,z+1).u[1];
      const Real dwdz2 = l(x,y,z-1).u[2]-2.0*l(x,y,z).u[2]+ l(x,y,z+1).u[2];

      Real dudxdy1, dvdxdy1, dwdxdy1;
      Real dudxdz1, dvdxdz1, dwdxdz1;
      Real dudydz1, dvdydz1, dwdydz1;
      if (inrange(x+2*sx) && inrange(y+2*sy))
      {
        dudxdy1 = sx*sy*(-0.5*( -1.5*l(x+2*sx,y,z).u[0]+2*l(x+2*sx,y+sy,z).u[0]-0.5*l(x+2*sx,y+2*sy,z).u[0] ) + 2*(-1.5*l(x+sx,y,z).u[0]+2*l(x+sx,y+sy,z).u[0]-0.5*l(x+sx,y+2*sy,z).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y+sy,z).u[0]-0.5*l(x,y+2*sy,z).u[0]));
        dvdxdy1 = sx*sy*(-0.5*( -1.5*l(x+2*sx,y,z).u[1]+2*l(x+2*sx,y+sy,z).u[1]-0.5*l(x+2*sx,y+2*sy,z).u[1] ) + 2*(-1.5*l(x+sx,y,z).u[1]+2*l(x+sx,y+sy,z).u[1]-0.5*l(x+sx,y+2*sy,z).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y+sy,z).u[1]-0.5*l(x,y+2*sy,z).u[1]));
        dwdxdy1 = sx*sy*(-0.5*( -1.5*l(x+2*sx,y,z).u[2]+2*l(x+2*sx,y+sy,z).u[2]-0.5*l(x+2*sx,y+2*sy,z).u[2] ) + 2*(-1.5*l(x+sx,y,z).u[2]+2*l(x+sx,y+sy,z).u[2]-0.5*l(x+sx,y+2*sy,z).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y+sy,z).u[2]-0.5*l(x,y+2*sy,z).u[2]));        
      }
      else
      {
        dudxdy1 = sx*sy*(l(x+sx,y+sy,z).u[0]-l(x+sx,y,z).u[0]) - (l(x,y+sy,z).u[0]-l(x,y,z).u[0]);
        dvdxdy1 = sx*sy*(l(x+sx,y+sy,z).u[1]-l(x+sx,y,z).u[1]) - (l(x,y+sy,z).u[1]-l(x,y,z).u[1]);
        dwdxdy1 = sx*sy*(l(x+sx,y+sy,z).u[2]-l(x+sx,y,z).u[2]) - (l(x,y+sy,z).u[2]-l(x,y,z).u[2]);
      }
      if (inrange(x+2*sx) && inrange(z+2*sz))
      {
        dudxdz1 = sx*sz*(-0.5*( -1.5*l(x+2*sx,y,z).u[0]+2*l(x+2*sx,y,z+sz).u[0]-0.5*l(x+2*sx,y,z+2*sz).u[0] ) + 2*(-1.5*l(x+sx,y,z).u[0]+2*l(x+sx,y,z+sz).u[0]-0.5*l(x+sx,y,z+2*sz).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+sz).u[0]-0.5*l(x,y,z+2*sz).u[0]));
        dvdxdz1 = sx*sz*(-0.5*( -1.5*l(x+2*sx,y,z).u[1]+2*l(x+2*sx,y,z+sz).u[1]-0.5*l(x+2*sx,y,z+2*sz).u[1] ) + 2*(-1.5*l(x+sx,y,z).u[1]+2*l(x+sx,y,z+sz).u[1]-0.5*l(x+sx,y,z+2*sz).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+sz).u[1]-0.5*l(x,y,z+2*sz).u[1]));
        dwdxdz1 = sx*sz*(-0.5*( -1.5*l(x+2*sx,y,z).u[2]+2*l(x+2*sx,y,z+sz).u[2]-0.5*l(x+2*sx,y,z+2*sz).u[2] ) + 2*(-1.5*l(x+sx,y,z).u[2]+2*l(x+sx,y,z+sz).u[2]-0.5*l(x+sx,y,z+2*sz).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+sz).u[2]-0.5*l(x,y,z+2*sz).u[2]));        
      }
      else
      {
        dudxdz1 = sx*sz*(l(x+sx,y,z+sz).u[0]-l(x+sx,y,z).u[0]) - (l(x,y,z+sz).u[0]-l(x,y,z).u[0]);
        dvdxdz1 = sx*sz*(l(x+sx,y,z+sz).u[1]-l(x+sx,y,z).u[1]) - (l(x,y,z+sz).u[1]-l(x,y,z).u[1]);
        dwdxdz1 = sx*sz*(l(x+sx,y,z+sz).u[2]-l(x+sx,y,z).u[2]) - (l(x,y,z+sz).u[2]-l(x,y,z).u[2]);
      }
      if (inrange(y+2*sy) && inrange(z+2*sz))
      {
        dudydz1 = sy*sz*(-0.5*( -1.5*l(x,y+2*sy,z).u[0]+2*l(x,y+2*sy,z+sz).u[0]-0.5*l(x,y+2*sy,z+2*sz).u[0] ) + 2*(-1.5*l(x,y+sy,z).u[0]+2*l(x,y+sy,z+sz).u[0]-0.5*l(x,y+sy,z+2*sz).u[0]) -1.5*(-1.5*l(x,y,z).u[0]+2*l(x,y,z+sz).u[0]-0.5*l(x,y,z+2*sz).u[0]));
        dvdydz1 = sy*sz*(-0.5*( -1.5*l(x,y+2*sy,z).u[1]+2*l(x,y+2*sy,z+sz).u[1]-0.5*l(x,y+2*sy,z+2*sz).u[1] ) + 2*(-1.5*l(x,y+sy,z).u[1]+2*l(x,y+sy,z+sz).u[1]-0.5*l(x,y+sy,z+2*sz).u[1]) -1.5*(-1.5*l(x,y,z).u[1]+2*l(x,y,z+sz).u[1]-0.5*l(x,y,z+2*sz).u[1]));
        dwdydz1 = sy*sz*(-0.5*( -1.5*l(x,y+2*sy,z).u[2]+2*l(x,y+2*sy,z+sz).u[2]-0.5*l(x,y+2*sy,z+2*sz).u[2] ) + 2*(-1.5*l(x,y+sy,z).u[2]+2*l(x,y+sy,z+sz).u[2]-0.5*l(x,y+sy,z+2*sz).u[2]) -1.5*(-1.5*l(x,y,z).u[2]+2*l(x,y,z+sz).u[2]-0.5*l(x,y,z+2*sz).u[2]));        
      }
      else
      {
        dudydz1 = sy*sz*(l(x,y+sy,z+sz).u[0]-l(x,y+sy,z).u[0]) - (l(x,y,z+sz).u[0]-l(x,y,z).u[0]);
        dvdydz1 = sy*sz*(l(x,y+sy,z+sz).u[1]-l(x,y+sy,z).u[1]) - (l(x,y,z+sz).u[1]-l(x,y,z).u[1]);
        dwdydz1 = sy*sz*(l(x,y+sy,z+sz).u[2]-l(x,y+sy,z).u[2]) - (l(x,y,z+sz).u[2]-l(x,y,z).u[2]);
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

  KernelComputeForces K(sim);
  cubism::compute<KernelComputeForces,VectorGrid,VectorLab,ScalarGrid,ScalarLab>(K,*sim.vel,*sim.chi);
  // do the final reductions and so on
  sim.obstacle_vector->computeForces();
}

CubismUP_3D_NAMESPACE_END