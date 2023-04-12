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

      const int sx = normX > 0 ? +1:-1;
      const int sy = normY > 0 ? +1:-1;
      const int sz = normZ > 0 ? +1:-1;

      VectorElement dveldx;
      if      (inrange(x+5*sx)) dveldx = sx*(  c0*l(x,y,z)+ c1*l(x+sx,y,z)+ c2*l(x+2*sx,y,z)+c3*l(x+3*sx,y,z)+c4*l(x+4*sx,y,z)+c5*l(x+5*sx,y,z));
      else if (inrange(x+2*sx)) dveldx = sx*(-1.5*l(x,y,z)+2.0*l(x+sx,y,z)-0.5*l(x+2*sx,y,z));
      else                      dveldx = sx*(l(x+sx,y,z)-l(x,y,z));
      VectorElement dveldy;
      if      (inrange(y+5*sy)) dveldy = sy*(  c0*l(x,y,z)+ c1*l(x,y+sy,z)+ c2*l(x,y+2*sy,z)+c3*l(x,y+3*sy,z)+c4*l(x,y+4*sy,z)+c5*l(x,y+5*sy,z));
      else if (inrange(y+2*sy)) dveldy = sy*(-1.5*l(x,y,z)+2.0*l(x,y+sy,z)-0.5*l(x,y+2*sy,z));
      else                      dveldy = sx*(l(x,y+sy,z)-l(x,y,z));
      VectorElement dveldz;
      if      (inrange(z+5*sz)) dveldz = sz*(  c0*l(x,y,z)+ c1*l(x,y,z+sz)+ c2*l(x,y,z+2*sz)+c3*l(x,y,z+3*sz)+c4*l(x,y,z+4*sz)+c5*l(x,y,z+5*sz));
      else if (inrange(z+2*sz)) dveldz = sz*(-1.5*l(x,y,z)+2.0*l(x,y,z+sz)-0.5*l(x,y,z+2*sz));
      else                      dveldz = sz*(l(x,y,z+sz)-l(x,y,z));

      const VectorElement dveldx2 = l(x-1,y,z)-2.0*l(x,y,z)+ l(x+1,y,z);
      const VectorElement dveldy2 = l(x,y-1,z)-2.0*l(x,y,z)+ l(x,y+1,z);
      const VectorElement dveldz2 = l(x,y,z-1)-2.0*l(x,y,z)+ l(x,y,z+1);

      VectorElement dveldxdy;
      VectorElement dveldxdz;
      VectorElement dveldydz;
      if (inrange(x+2*sx) && inrange(y+2*sy)) dveldxdy = sx*sy*(-0.5*( -1.5*l(x+2*sx,y     ,z     )+2*l(x+2*sx,y+  sy,z     )  -0.5*l(x+2*sx,y+2*sy,z     )         ) + 2*(-1.5*l(x+sx,y,z)+2*l(x+sx,y+sy,z)-0.5*l(x+sx,y+2*sy,z)) -1.5*(-1.5*l(x,y,z)+2*l(x,y+sy,z)-0.5*l(x,y+2*sy,z)));
      else                                    dveldxdy = sx*sy*(            l(x+  sx,y+  sy,z     )-  l(x+  sx,y     ,z     )) -   (l(x     ,y  +sy,z     )-l(x,y,z));
      if (inrange(y+2*sy) && inrange(z+2*sz)) dveldydz = sy*sz*(-0.5*( -1.5*l(x     ,y+2*sy,z     )+2*l(x     ,y+2*sy,z+  sz)  -0.5*l(x     ,y+2*sy,z+2*sz)         ) + 2*(-1.5*l(x,y+sy,z)+2*l(x,y+sy,z+sz)-0.5*l(x,y+sy,z+2*sz)) -1.5*(-1.5*l(x,y,z)+2*l(x,y,z+sz)-0.5*l(x,y,z+2*sz)));
      else                                    dveldydz = sy*sz*(            l(x     ,y+  sy,z+  sz)-  l(x     ,y+  sy,z     )) -   (l(x     ,y     ,z+  sz)-l(x,y,z));
      if (inrange(x+2*sx) && inrange(z+2*sz)) dveldxdz = sx*sz*(-0.5*( -1.5*l(x     ,y     ,z+2*sz)+2*l(x+  sx,y     ,z+2*sz)  -0.5*l(x+2*sx,y     ,z+2*sz)         ) + 2*(-1.5*l(x,y,z+sz)+2*l(x+sx,y,z+sz)-0.5*l(x+2*sx,y,z+sz)) -1.5*(-1.5*l(x,y,z)+2*l(x+sx,y,z)-0.5*l(x+2*sx,y,z)));
      else                                    dveldxdz = sx*sz*(            l(x+  sx,y     ,z+  sz)-  l(x     ,y     ,z+  sz)) -   (l(x  +sx,y     ,z     )-l(x,y,z));

      const Real dudx = dveldx.u[0] + dveldx2.u[0]*(ix-x) + dveldxdy.u[0]*(iy-y) + dveldxdz.u[0]*(iz-z);
      const Real dvdx = dveldx.u[1] + dveldx2.u[1]*(ix-x) + dveldxdy.u[1]*(iy-y) + dveldxdz.u[1]*(iz-z);
      const Real dwdx = dveldx.u[2] + dveldx2.u[2]*(ix-x) + dveldxdy.u[2]*(iy-y) + dveldxdz.u[2]*(iz-z);
      const Real dudy = dveldy.u[0] + dveldy2.u[0]*(iy-y) + dveldydz.u[0]*(iz-z) + dveldxdy.u[0]*(ix-x);
      const Real dvdy = dveldy.u[1] + dveldy2.u[1]*(iy-y) + dveldydz.u[1]*(iz-z) + dveldxdy.u[1]*(ix-x);
      const Real dwdy = dveldy.u[2] + dveldy2.u[2]*(iy-y) + dveldydz.u[2]*(iz-z) + dveldxdy.u[2]*(ix-x);
      const Real dudz = dveldz.u[0] + dveldz2.u[0]*(iz-z) + dveldxdz.u[0]*(ix-x) + dveldydz.u[0]*(iy-y);
      const Real dvdz = dveldz.u[1] + dveldz2.u[1]*(iz-z) + dveldxdz.u[1]*(ix-x) + dveldydz.u[1]*(iy-y);
      const Real dwdz = dveldz.u[2] + dveldz2.u[2]*(iz-z) + dveldxdz.u[2]*(ix-x) + dveldydz.u[2]*(iy-y);

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
      o->P[i] = P;
      o->fX[i] = -P * dx + _1oH * (dudx * dx + dudy * dy + dudz * dz);
      o->fY[i] = -P * dy + _1oH * (dvdx * dx + dvdy * dy + dvdz * dz);
      o->fZ[i] = -P * dz + _1oH * (dwdx * dx + dwdy * dy + dwdz * dz);
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
