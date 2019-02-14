//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorAdvectDiffuse_h
#define CubismUP_3D_CoordinatorAdvectDiffuse_h

#include "PenalizationObstacleVisitor.h"

class OperatorMinusDivTmpU : public GenericLabOperator
{
 private:
  const double dt;
 public:
  OperatorMinusDivTmpU(double _dt) : dt(_dt)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 4, 0,4,5,6,7);
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
  }
  ~OperatorMinusDivTmpU() {}
  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = 0.5 * info.h_gridpoint*info.h_gridpoint / dt, DT = dt;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // Poisson solver reads field p for the rhs
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real dXx_dux = (LE.chi-LW.chi)*( L.tmpU + DT*(LE.p-LW.p) );
      const Real dXy_duy = (LN.chi-LS.chi)*( L.tmpV + DT*(LN.p-LS.p) );
      const Real dXz_duz = (LB.chi-LF.chi)*( L.tmpW + DT*(LB.p-LF.p) );
      const Real divU = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
      const Real dU_XX = L.chi * L.chi * divU, inv1pX = 1 / (1 + L.chi);
      o(ix,iy,iz).p = - fac*dU_XX + fac*inv1pX*(dXx_dux+dXy_duy+dXz_duz);
    }
  }
};

class OperatorAdvectDiffuse : public GenericLabOperator
{
  private:
  const double dt;
  const double mu;
  const double lambda;
  const Real* const uInf;

  public:
   OperatorAdvectDiffuse(const double _dt, double m, const double l,
     const Real* const u) : dt(_dt), mu(m), lambda(l), uInf(u)
   {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2; stencil_end[1] = 2;  stencil_end[2] = 2;
   }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real facA = -dt/(2*info.h_gridpoint);
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real dudx= LE.u-LW.u, dvdx= LE.v-LW.v, dwdx= LE.w-LW.w;
      const Real dudy= LN.u-LS.u, dvdy= LN.v-LS.v, dwdy= LN.w-LS.w;
      const Real dudz= LB.u-LF.u, dvdz= LB.v-LF.v, dwdz= LB.w-LF.w;
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
      const Real duD = LN.u+LS.u + LE.u+LW.u + LF.u+LB.u - L.u*6;
      const Real dvD = LN.v+LS.v + LE.v+LW.v + LF.v+LB.v - L.v*6;
      const Real dwD = LN.w+LS.w + LE.w+LW.w + LF.w+LB.w - L.w*6;
      const Real duA = u * dudx + v * dudy + w * dudz;
      const Real dvA = u * dvdx + v * dvdy + w * dvdz;
      const Real dwA = u * dwdx + v * dwdy + w * dwdz;
      o(ix,iy,iz).tmpU = L.u + facA*duA + facD*duD;
      o(ix,iy,iz).tmpV = L.v + facA*dvA + facD*dvD;
      o(ix,iy,iz).tmpW = L.w + facA*dwA + facD*dwD;
    }
  }
};

template <typename Lab>
class CoordinatorAdvectDiffuse : public GenericCoordinator
{
public:
  CoordinatorAdvectDiffuse(SimulationData & s) : GenericCoordinator(s) { }

  ~CoordinatorAdvectDiffuse() { }

  void operator()(const double dt)
  {
    check("AdvectDiffuse - start");
    const int nthreads = omp_get_max_threads();
    {
      //zero fields, going to contain Udef:
      #pragma omp parallel for schedule(static)
      for(unsigned i=0; i<vInfo.size(); i++) {
        const BlockInfo& info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
          b(ix,iy,iz).tmpU = 0; b(ix,iy,iz).tmpV = 0; b(ix,iy,iz).tmpW = 0;
        }
      }
      //store deformation velocities onto tmp fields:
      ObstacleVisitor*visitor=new PenalizationObstacleVisitor(grid,dt,sim.uinf);
      sim.obstacle_vector->Accept(visitor);
      delete visitor;
    }
    {   //place onto p: ( div u^(t+1) - div u^* ) / dt
      //where i want div u^(t+1) to be equal to div udef
      std::vector<OperatorMinusDivTmpU*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i) diff[i] = new OperatorMinusDivTmpU(dt);

      compute<OperatorMinusDivTmpU>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }
    {
      using advection = OperatorAdvectDiffuse;
      std::vector<advection*> adv1(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0; i<nthreads; ++i)
        adv1[i] = new advection(dt, sim.nu, sim.lambda, sim.uinf);
      compute(adv1);
      for(int i=0; i<nthreads; i++) delete adv1[i];
    }

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      const BlockInfo& info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
        b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
        b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
      }
    }

    check("AdvectDiffuse - end");
  }

  std::string getName()
  {
    return "AdvectDiffuse";
  }
};
#endif
