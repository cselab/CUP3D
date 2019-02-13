//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorAdvectDiffuse_h
#define CubismUP_3D_CoordinatorAdvectDiffuse_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"
#include <cmath>

struct PenalizationObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const double dt;
  const Real * const uInf;
  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();

  PenalizationObstacleVisitor(FluidGridMPI*g, const double _dt,
    const Real*const u) : grid(g), dt(_dt), uInf(u) { }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    //using CHI_MAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
    using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
    #pragma omp parallel
    {
      const std::map<int, ObstacleBlock*>& obstblocks =
        obstacle->getObstacleBlocks();
      double uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);
      const size_t Nblocks = vInfo.size();
      #pragma omp for schedule(dynamic)
      for (size_t i = 0; i < Nblocks; ++i)
      {
        const BlockInfo& info = vInfo[i];
        const auto pos = obstblocks.find(info.blockID);
        if(pos == obstblocks.end()) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        UDEFMAT & __restrict__ UDEF = pos->second->udef;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          Real p[3]; info.pos(p, ix, iy, iz);
          p[0]-=centerOfMass[0]; p[1]-=centerOfMass[1]; p[2]-=centerOfMass[2];
          const Real object_UR[3] = {
              (Real) omegaBody[1]*p[2] - (Real) omegaBody[2]*p[1],
              (Real) omegaBody[2]*p[0] - (Real) omegaBody[0]*p[2],
              (Real) omegaBody[0]*p[1] - (Real) omegaBody[1]*p[0]
          };
          const Real U_TOT[3] = {
              (Real)uBody[0] +object_UR[0] +UDEF[iz][iy][ix][0] -uInf[0],
              (Real)uBody[1] +object_UR[1] +UDEF[iz][iy][ix][1] -uInf[1],
              (Real)uBody[2] +object_UR[2] +UDEF[iz][iy][ix][2] -uInf[2]
          };
          // what if multiple obstacles share a block??
          // let's plus equal and wake up during the night to stress about it
          b(ix,iy,iz).tmpU += U_TOT[0];
          b(ix,iy,iz).tmpV += U_TOT[1];
          b(ix,iy,iz).tmpW += U_TOT[2];
        }
      }
    }
  }
};

class OperatorMinusDivTmpU : public GenericLabOperator
{
 private:
  const double dt;
 public:
  OperatorMinusDivTmpU(double _dt) : dt(_dt)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 5,6,7);
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
  }
  ~OperatorMinusDivTmpU() {}
  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = 0.5 * info.h_gridpoint*info.h_gridpoint / dt;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      // Poisson solver reads field p for the rhs
      const Real dUdef = lab(ix+1,iy  ,iz  ).tmpU - lab(ix-1,iy  ,iz  ).tmpU;
      const Real dVdef = lab(ix  ,iy+1,iz  ).tmpV - lab(ix  ,iy-1,iz  ).tmpV;
      const Real dWdef = lab(ix  ,iy  ,iz+1).tmpW - lab(ix  ,iy  ,iz-1).tmpW;
      o(ix,iy,iz).p = - fac * lab(ix,iy,iz).chi * (dUdef+dVdef+dWdef);
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
    const Real facP = dt * lambda;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);

      const Real US = lab(ix,iy,iz).tmpU, VS = lab(ix,iy,iz).tmpV;
      const Real WS = lab(ix,iy,iz).tmpW, CHI = lab(ix,iy,iz).chi;
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
      const Real duP = facP*CHI * (US - L.u) / (1 + CHI*facP);
      const Real dvP = facP*CHI * (VS - L.v) / (1 + CHI*facP);
      const Real dwP = facP*CHI * (WS - L.w) / (1 + CHI*facP);
      o(ix,iy,iz).tmpU = L.u + facA*duA + facD*duD + duP;
      o(ix,iy,iz).tmpV = L.v + facA*dvA + facD*dvD + dvP;
      o(ix,iy,iz).tmpW = L.w + facA*dwA + facD*dwD + dwP;
    }
  }
};


template <typename Lab>
class CoordinatorAdvectDiffuse : public GenericCoordinator
{
protected:
  const Real MU;
  double* const lambda;
  const Real* const uInf;
  IF3D_ObstacleVector** const obstacleVector;

public:
  CoordinatorAdvectDiffuse(const Real m, double*const l, const Real*const u,
    IF3D_ObstacleVector** const ov, FluidGridMPI * g) : GenericCoordinator(g),
    MU(m), lambda(l), uInf(u), obstacleVector(ov) { }

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
      ObstacleVisitor*visitor = new PenalizationObstacleVisitor(grid,dt,uInf);
      (*obstacleVector)->Accept(visitor);
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
      for(int i=0;i<nthreads;++i) adv1[i] = new advection(dt,MU,*lambda,uInf);

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
