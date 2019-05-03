//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/IterativePressurePenalization.h"
#include "obstacles/ObstacleVector.h"
#ifdef _ACCFFT_
#include "poisson/PoissonSolverACCPeriodic.h"
#include "poisson/PoissonSolverACCUnbounded.h"
#else
#include "poisson/PoissonSolverPeriodic.h"
#include "poisson/PoissonSolverUnbounded.h"
#endif
// TODO : Cosine transform on GPU!?
#include "poisson/PoissonSolverMixed.h"
#include "poisson/PoissonSolverHYPREMixed.h"
#include "poisson/PoissonSolverPETSCMixed.h"

// define this to update obstacles with old (mrag-like) approach of integrating
// momenta contained in chi before the penalization step:
#define EXPL_INTEGRATE_MOM

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
static constexpr Real DBLEPS = std::numeric_limits<double>::epsilon();

namespace
{

static inline PenalizationBlock* getPenalBlockPtr(
  PenalizationGridMPI*const grid, const int blockID) {
  assert(grid not_eq nullptr);
  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();
  return (PenalizationBlock*) vInfo[blockID].ptrBlock;
}

struct PressureRHSObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();

  PressureRHSObstacleVisitor(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obstacle)
  {
    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        const cubism::BlockInfo& info = vInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        const UDEFMAT & __restrict__ UDEF = pos->udef;
        const  CHIMAT & __restrict__  CHI = pos->chi;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          // What if multiple obstacles share a block? Do not write udef onto
          // grid if CHI stored on the grid is greater than obst's CHI.
          if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
          // What if two obstacles overlap? Let's plus equal. After all here
          // we are computing divUs, maybe one obstacle has divUs 0. We will
          // need a repulsion term of the velocity at some point in the code.
          b(ix,iy,iz).tmpU += UDEF[iz][iy][ix][0];
          b(ix,iy,iz).tmpV += UDEF[iz][iy][ix][1];
          b(ix,iy,iz).tmpW += UDEF[iz][iy][ix][2];
        }
      }
    }
  }
};

struct KernelPressureRHS
{
  Real dt, lamdt;
  const Real fadeLen[3], ext[3], iFade[3];
  static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
  PoissonSolver * const solver;

  inline bool _is_touching(const FluidBlock& b) const {
    const bool touchW = fadeLen[0] >= b.min_pos[0];
    const bool touchE = fadeLen[0] >= ext[0] - b.max_pos[0];
    const bool touchS = fadeLen[1] >= b.min_pos[1];
    const bool touchN = fadeLen[1] >= ext[1] - b.max_pos[1];
    const bool touchB = fadeLen[2] >= b.min_pos[2];
    const bool touchF = fadeLen[2] >= ext[2] - b.max_pos[2];
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }

  inline Real fade(const BlockInfo&i, const int x,const int y,const int z) const
  {
    Real p[3]; i.pos(p, x, y, z);
    const Real zt = iFade[2] * std::max(Real(0), fadeLen[2] -(ext[2]-p[2]) );
    const Real zb = iFade[2] * std::max(Real(0), fadeLen[2] - p[2] );
    const Real yt = iFade[1] * std::max(Real(0), fadeLen[1] -(ext[1]-p[1]) );
    const Real yb = iFade[1] * std::max(Real(0), fadeLen[1] - p[1] );
    const Real xt = iFade[0] * std::max(Real(0), fadeLen[0] -(ext[0]-p[0]) );
    const Real xb = iFade[0] * std::max(Real(0), fadeLen[0] - p[0] );
    return 1-std::pow(std::min( std::max({zt,zb,yt,yb,xt,xb}), (Real)1), 2);
  }

  inline Real RHS(Lab&l, const int x,const int y,const int z) const
  {
    const FluidElement & L  = l(x,  y,  z);
    const FluidElement & LW = l(x-1,y,  z  ), & LE = l(x+1,y,  z  );
    const FluidElement & LS = l(x,  y-1,z  ), & LN = l(x,  y+1,z  );
    const FluidElement & LF = l(x,  y,  z-1), & LB = l(x,  y,  z+1);
    const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
    const Real divUf = LE.u-LW.u + LN.v-LS.v + LB.w-LF.w;
    return divUf + L.chi * divUs;
  }

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1,2,2,2,false,6,1,2,3,5,6,7);

  KernelPressureRHS(double _dt, double L, const Real B[3], std::array<Real,3> E,
    PoissonSolver* ps) : dt(_dt), lamdt(_dt*L), fadeLen{B[0],B[1],B[2]},
    ext{E[0],E[1],E[2]}, iFade{1/(B[0]+EPS),1/(B[1]+EPS),1/(B[2]+EPS)},
    solver(ps) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt;
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;
    if( not _is_touching(o) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        ret[SZ*iz + SY*iy + SX*ix] = fac * RHS(lab, ix,iy,iz);
        //o(ix,iy,iz).p = ret[SZ*iz + SY*iy + SX*ix];
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        ret[SZ*iz + SY*iy + SX*ix] = fac*fade(info,ix,iy,iz)*RHS(lab,ix,iy,iz);
        //o(ix,iy,iz).p = ret[SZ*iz + SY*iy + SX*ix];
      }
    }
  }
};

// place in PenalizationBlocks velocity before iterations plus pressure correct
struct KernelIterateGradP
{
  const Real dt;
  PenalizationGridMPI * const penGrid;
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);

  KernelIterateGradP(double _dt, PenalizationGridMPI * const pen) :
    dt(_dt), penGrid{pen} {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
       // p contains the pressure correction after the Poisson solver
      const Real dU = fac * ( lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p );
      const Real dV = fac * ( lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p );
      const Real dW = fac * ( lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p );
      t(ix,iy,iz).uPres = lab(ix,iy,iz).u + dU;
      t(ix,iy,iz).vPres = lab(ix,iy,iz).v + dV;
      t(ix,iy,iz).wPres = lab(ix,iy,iz).w + dW;
    }
  }
};

// place on grid the penalized velocity plus pressure correction
struct KernelGradP
{
  const Real dt;
  PenalizationGridMPI * const penGrid;
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);

  KernelGradP(double _dt, PenalizationGridMPI * const pen): dt(_dt), penGrid{pen} {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    const PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
       const Real dU = fac * ( lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p );
       const Real dV = fac * ( lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p );
       const Real dW = fac * ( lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p );
       o(ix,iy,iz).u = t(ix,iy,iz).uPenl + dU;
       o(ix,iy,iz).v = t(ix,iy,iz).vPenl + dV;
       o(ix,iy,iz).w = t(ix,iy,iz).wPenl + dW;
    }
  }
};

struct KernelIntegrateFluidMomenta : public ObstacleVisitor
{
  const double lambda, dt;
  ObstacleVector * const obstacle_vector;
  PenalizationGridMPI * const penGrid;
  const cubism::BlockInfo * info_ptr = nullptr;
  inline double dvol(const cubism::BlockInfo&info, const int x, const int y, const int z) const {
    double h[3]; info.spacing(h, x, y, z);
    return h[0] * h[1] * h[2];
  }

  KernelIntegrateFluidMomenta(double _dt, double _lambda, ObstacleVector* ov,
    PenalizationGridMPI*const pen) : lambda(_lambda), dt(_dt),
    obstacle_vector(ov), penGrid(pen) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr);
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
  }

  void visit(Obstacle* const op)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const std::vector<ObstacleBlock*>& obstblocks = op->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const std::array<double,3> CM = op->getCenterOfMass();
    const PenalizationBlock& b = * getPenalBlockPtr(penGrid, info.blockID);

    const CHIMAT & __restrict__ CHI = o->chi;
    double &VV = o->V;
    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    double &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    double &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;

    #ifndef EXPL_INTEGRATE_MOM
      const UDEFMAT & __restrict__ UDEF = o->udef;
      o->GfX = 0;
      o->GpX = 0; o->GpY = 0; o->GpZ = 0;
      o->Gj0 = 0; o->Gj1 = 0; o->Gj2 = 0;
      o->Gj3 = 0; o->Gj4 = 0; o->Gj5 = 0;
      o->GuX = 0; o->GuY = 0; o->GuZ = 0;
      o->GaX = 0; o->GaY = 0; o->GaZ = 0;
      const Real lambdt = lambda*dt;
    #endif

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (CHI[iz][iy][ix] <= 0) continue;
      double p[3]; info.pos(p, ix, iy, iz);
      const double dv = dvol(info, ix, iy, iz), X = CHI[iz][iy][ix];
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      VV += X * dv;
      J0 += X * dv * ( p[1]*p[1] + p[2]*p[2] );
      J1 += X * dv * ( p[0]*p[0] + p[2]*p[2] );
      J2 += X * dv * ( p[0]*p[0] + p[1]*p[1] );
      J3 -= X * dv * p[0]*p[1];
      J4 -= X * dv * p[0]*p[2];
      J5 -= X * dv * p[1]*p[2];

      FX += X * dv * b(ix,iy,iz).uPres;
      FY += X * dv * b(ix,iy,iz).vPres;
      FZ += X * dv * b(ix,iy,iz).wPres;
      TX += X * dv * ( p[1] * b(ix,iy,iz).wPres - p[2] * b(ix,iy,iz).vPres );
      TY += X * dv * ( p[2] * b(ix,iy,iz).uPres - p[0] * b(ix,iy,iz).wPres );
      TZ += X * dv * ( p[0] * b(ix,iy,iz).vPres - p[1] * b(ix,iy,iz).uPres );

      #ifndef EXPL_INTEGRATE_MOM
        const Real penalFac = dv * lambdt * X / ( 1 + X * lambdt );
        o->GfX += penalFac;
        o->GpX += penalFac * p[0];
        o->GpY += penalFac * p[1];
        o->GpZ += penalFac * p[2];
        o->Gj0 += penalFac * ( p[1]*p[1] + p[2]*p[2] );
        o->Gj1 += penalFac * ( p[0]*p[0] + p[2]*p[2] );
        o->Gj2 += penalFac * ( p[0]*p[0] + p[1]*p[1] );
        o->Gj3 -= penalFac * p[0]*p[1];
        o->Gj4 -= penalFac * p[0]*p[2];
        o->Gj5 -= penalFac * p[1]*p[2];
        const double DiffU[3] = {
          b(ix,iy,iz).uPres - UDEF[iz][iy][ix][0],
          b(ix,iy,iz).vPres - UDEF[iz][iy][ix][1],
          b(ix,iy,iz).wPres - UDEF[iz][iy][ix][2]
        };
        o->GuX += penalFac * DiffU[0];
        o->GuY += penalFac * DiffU[1];
        o->GuZ += penalFac * DiffU[2];
        o->GaX += penalFac * ( p[1] * DiffU[2] - p[2] * DiffU[1] );
        o->GaY += penalFac * ( p[2] * DiffU[0] - p[0] * DiffU[2] );
        o->GaZ += penalFac * ( p[0] * DiffU[1] - p[1] * DiffU[0] );
      #endif
    }
  }
};

struct KernelFinalizeObstacleVel : public ObstacleVisitor
{
  const double dt, lambda;
  FluidGridMPI * const grid;

  KernelFinalizeObstacleVel(double _dt, double _lambda, FluidGridMPI*g) :
    dt(_dt), lambda(_lambda), grid(g) { }

  void visit(Obstacle* const obst)
  {
    static constexpr int nQoI = 29;
    double M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); i++)
    {
      if(oBlock[i] == nullptr) continue;
      int k = 0;
      M[k++] += oBlock[i]->V ;
      M[k++] += oBlock[i]->FX; M[k++] += oBlock[i]->FY; M[k++] += oBlock[i]->FZ;
      M[k++] += oBlock[i]->TX; M[k++] += oBlock[i]->TY; M[k++] += oBlock[i]->TZ;
      M[k++] += oBlock[i]->J0; M[k++] += oBlock[i]->J1; M[k++] += oBlock[i]->J2;
      M[k++] += oBlock[i]->J3; M[k++] += oBlock[i]->J4; M[k++] += oBlock[i]->J5;
      #ifndef EXPL_INTEGRATE_MOM
        M[k++]+=oBlock[i]->GfX;
        M[k++]+=oBlock[i]->GpX; M[k++]+=oBlock[i]->GpY; M[k++]+=oBlock[i]->GpZ;
        M[k++]+=oBlock[i]->Gj0; M[k++]+=oBlock[i]->Gj1; M[k++]+=oBlock[i]->Gj2;
        M[k++]+=oBlock[i]->Gj3; M[k++]+=oBlock[i]->Gj4; M[k++]+=oBlock[i]->Gj5;
        M[k++]+=oBlock[i]->GuX; M[k++]+=oBlock[i]->GuY; M[k++]+=oBlock[i]->GuZ;
        M[k++]+=oBlock[i]->GaX; M[k++]+=oBlock[i]->GaY; M[k++]+=oBlock[i]->GaZ;
        assert(k==29);
      #else
        assert(k==13);
      #endif
    }
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_DOUBLE, MPI_SUM, comm);
    assert(std::fabs(obst->mass - M[ 0]) < 10*DBLEPS);
    assert(std::fabs(obst->J[0] - M[ 7]) < 10*DBLEPS);
    assert(std::fabs(obst->J[1] - M[ 8]) < 10*DBLEPS);
    assert(std::fabs(obst->J[2] - M[ 9]) < 10*DBLEPS);
    assert(std::fabs(obst->J[3] - M[10]) < 10*DBLEPS);
    assert(std::fabs(obst->J[4] - M[11]) < 10*DBLEPS);
    assert(std::fabs(obst->J[5] - M[12]) < 10*DBLEPS);
    assert(M[0] > DBLEPS);

    #ifndef EXPL_INTEGRATE_MOM
      obst->penalM    = M[13];
      obst->penalCM   = { M[14], M[15], M[16] };
      obst->penalJ    = { M[17], M[18], M[19], M[20], M[21], M[22] };
      obst->penalLmom = { M[23], M[24], M[25] };
      obst->penalAmom = { M[26], M[27], M[28] };
    #else
      obst->penalM    = M[0];
      obst->penalCM   = { 0, 0, 0 };
      obst->penalJ    = { M[ 7], M[ 8], M[ 9], M[10], M[11], M[12] };
      obst->penalLmom = { M[1], M[2], M[3] };
      obst->penalAmom = { M[4], M[5], M[6] };
    #endif

    obst->computeVelocities();
  }
};

struct KernelFinalizePenalizationForce : public ObstacleVisitor
{
  FluidGridMPI * const grid;

  KernelFinalizePenalizationForce(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obst)
  {
    static constexpr int nQoI = 6;
    double M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); ++i) {
      if(oBlock[i] == nullptr) continue;
      M[0] += oBlock[i]->FX; M[1] += oBlock[i]->FY; M[2] += oBlock[i]->FZ;
      M[3] += oBlock[i]->TX; M[4] += oBlock[i]->TY; M[5] += oBlock[i]->TZ;
    }
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_DOUBLE, MPI_SUM, comm);
    obst->force[0]  = M[0];
    obst->force[1]  = M[1];
    obst->force[2]  = M[2];
    obst->torque[0] = M[3];
    obst->torque[1] = M[4];
    obst->torque[2] = M[5];
  }
};

struct KernelPenalization : public ObstacleVisitor
{
  Real MX = 0, MY = 0, MZ = 0, DMX = 0, DMY = 0, DMZ = 0;
  const double lambda, dt; const int iter;
  ObstacleVector * const obstacle_vector;
  PenalizationGridMPI * const penGrid;
  PenalizationGridMPI * const accGrid;
  const cubism::BlockInfo * info_ptr = nullptr;

  KernelPenalization(double _lambda, double _dt, ObstacleVector* ov, int _iter,
    PenalizationGridMPI* const _pen, PenalizationGridMPI* const _acc) :
    lambda(_lambda), dt(_dt), iter(_iter), obstacle_vector(ov), penGrid(_pen), accGrid(_acc) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);

    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;

    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();
    const double dv = std::pow(info.h_gridpoint, 3);

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // What if multiple obstacles share a block? Do not write udef onto
      // grid if CHI stored on the grid is greater than obst's CHI.
      if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      if(CHI[iz][iy][ix] <= 0) continue; // no need to do anything
      double p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];
      const double U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };
      #ifdef EXPL_INTEGRATE_MOM
        const Real penalFac = lambda;
      #else
        const Real penalFac = lambda / (1 + dt * CHI[iz][iy][ix] * lambda);
      #endif
      const Real FPX = penalFac * (U_TOT[0] - t(ix,iy,iz).uPres);
      const Real FPY = penalFac * (U_TOT[1] - t(ix,iy,iz).vPres);
      const Real FPZ = penalFac * (U_TOT[2] - t(ix,iy,iz).wPres);

      const Real UPenal = b(ix,iy,iz).u + dt * CHI[iz][iy][ix] * FPX;
      const Real VPenal = b(ix,iy,iz).v + dt * CHI[iz][iy][ix] * FPY;
      const Real WPenal = b(ix,iy,iz).w + dt * CHI[iz][iy][ix] * FPZ;
      const Real DPX = UPenal - t(ix,iy,iz).uPenl;
      const Real DPY = VPenal - t(ix,iy,iz).vPenl;
      const Real DPZ = WPenal - t(ix,iy,iz).wPenl;

      t(ix,iy,iz).uPenl = UPenal;
      t(ix,iy,iz).vPenl = VPenal;
      t(ix,iy,iz).wPenl = WPenal;
      b(ix,iy,iz).tmpU = dt*FPX - UDEF[iz][iy][ix][0]; // used to compute press
      b(ix,iy,iz).tmpV = dt*FPY - UDEF[iz][iy][ix][1]; // rhs: includes penal
      b(ix,iy,iz).tmpW = dt*FPZ - UDEF[iz][iy][ix][2]; // force-udef, without X

      FX += dv * CHI[iz][iy][ix] * FPX;
      FY += dv * CHI[iz][iy][ix] * FPY;
      FZ += dv * CHI[iz][iy][ix] * FPZ;
      TX += dv * CHI[iz][iy][ix] * ( p[1] * FPZ - p[2] * FPY );
      TY += dv * CHI[iz][iy][ix] * ( p[2] * FPX - p[0] * FPZ );
      TZ += dv * CHI[iz][iy][ix] * ( p[0] * FPY - p[1] * FPX );
      MX += std::pow(UPenal, 2); DMX += std::pow(DPX, 2);
      MY += std::pow(VPenal, 2); DMY += std::pow(DPY, 2);
      MZ += std::pow(WPenal, 2); DMZ += std::pow(DPZ, 2);
    }
  }
};

}

void IterativePressurePenalization::initializeFields()
{
  /*
  sim.startProfiler("PresRHS Udef");
  if(sim.obstacle_vector->nObstacles() > 0)
  { //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      assert((size_t) vInfo[i].blockID == i);
      FluidBlock&        b = *(FluidBlock*) vInfo[i].ptrBlock;
      PenalizationBlock& t = * getPenalBlockPtr(penalizationGrid, i);
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        t(ix,iy,iz).uPre = b(ix,iy,iz).u;
        t(ix,iy,iz).vPre = b(ix,iy,iz).v;
        t(ix,iy,iz).wPre = b(ix,iy,iz).w;
        b(ix,iy,iz).tmpU = 0;
        b(ix,iy,iz).tmpV = 0;
        b(ix,iy,iz).tmpW = 0;
      }
    }
    //store deformation velocities onto tmp fields:
    ObstacleVisitor* visitor = new PressureRHSObstacleVisitor(grid);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();
  */
}

IterativePressurePenalization::IterativePressurePenalization(SimulationData & s) : Operator(s)
{
  if(sim.bUseFourierBC)
  pressureSolver = new PoissonSolverPeriodic(sim);
  else if (sim.bUseUnboundedBC)
  pressureSolver = new PoissonSolverUnbounded(sim);
  #ifdef CUP_HYPRE
  else if (sim.useSolver == "hypre")
  pressureSolver = new PoissonSolverMixed_HYPRE(sim);
  #endif
  #ifdef CUP_PETSC
  else if (sim.useSolver == "petsc")
  pressureSolver = new PoissonSolverMixed_PETSC(sim);
  #endif
  else
  pressureSolver = new PoissonSolverMixed(sim);
  sim.pressureSolver = pressureSolver;

  // no need to allocate stretched mesh stuff here because we will never read
  // grid spacing from this grid!
  if(sim.rank==0) printf("Allocating the penalization helper grid.\n");

  penalizationGrid = new PenalizationGridMPI(
    sim.nprocsx, sim.nprocsy, sim.nprocsz, sim.local_bpdx,
    sim.local_bpdy, sim.local_bpdz, sim.maxextent, sim.app_comm);
}

void IterativePressurePenalization::operator()(const double dt)
{
  #ifdef EXPL_INTEGRATE_MOM
    sim.lambda = 1.0/dt;
  #endif
  // first copy velocity before either Pres or Penal onto penalization blocks
  // also put udef into tmpU fields
  initializeFields();

  int iter=0;
  Real relDF = 1e3;
  bool bConverged = false;
  for(iter = 0; iter < 1000; ++iter)
  {
    {
      // compute velocity after pressure projection (PP) but without Penal
      // using pressure computed at prior step / iteration
      sim.startProfiler("GradP"); //pressure correction dudt* = - grad P / rho
      const KernelIterateGradP K(dt, penalizationGrid);
      compute<KernelIterateGradP>(K);
      sim.stopProfiler();
    }

    sim.startProfiler("Obst Int Vel");
    { // integrate momenta by looping over grid
      #pragma omp parallel
      { // each thread needs to call its own non-const operator() function
        KernelIntegrateFluidMomenta K(dt, sim.lambda,
            sim.obstacle_vector, penalizationGrid);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
      }
    }
    sim.stopProfiler();

    sim.startProfiler("Obst Upd Vel");
    {
      ObstacleVisitor*K = new KernelFinalizeObstacleVel(dt, sim.lambda, sim.grid);
      sim.obstacle_vector->Accept(K); // accept you son of a french cow
      delete K;
    }
    sim.stopProfiler();

    // finally update vel with penalization but without pressure
    sim.startProfiler("Penalization");
    {
      double M[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      #pragma omp parallel reduction (+ : M[:6])
      { // each thread needs to call its own non-const operator() function
        KernelPenalization K(sim.lambda, dt, sim.obstacle_vector, iter, penalizationGrid, accelerationGrid);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
        M[0] += K.MX; M[3] += K.DMX;
        M[1] += K.MY; M[4] += K.DMY;
        M[2] += K.MZ; M[5] += K.DMZ;
      }
      const auto comm = sim.grid->getCartComm();
      MPI_Allreduce(MPI_IN_PLACE, M, 6, MPI_DOUBLE, MPI_SUM, comm);
      relDF = std::sqrt( (M[3]+M[4]+M[5]) / (EPS + M[0]+M[1]+M[2]) );

      ObstacleVisitor*K = new KernelFinalizePenalizationForce(sim.grid);
      sim.obstacle_vector->Accept(K); // accept you son of a french cow
      delete K;

    }
    sim.stopProfiler();

    if(sim.verbose) printf("iter:%02d - max relative error: %f\n", iter, relDF);
    //if(relDF < 0.001) break;
    if(iter>0 && relDF<0.001) bConverged = true;

    {
      sim.startProfiler("PresRHS Kernel");
      sim.pressureSolver->reset();
      //place onto p: ( div u^(t+1) - div u^* ) / dt
      //where i want div u^(t+1) to be equal to div udef
      const KernelPressureRHS K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
      compute<KernelPressureRHS>(K);
      sim.stopProfiler();
    }

    pressureSolver->solve();

    sim.startProfiler("sol2cub");
    pressureSolver->_fftw2cub();
    sim.stopProfiler();

    if(bConverged) break;
  }

  sim.startProfiler("GradP"); //pressure correction dudt* = - grad P / rho
  {
    const KernelGradP K(dt, penalizationGrid);
    compute<KernelGradP>(K);
  }
  sim.stopProfiler();

  /*
  if(not sim.muteAll)
  {
  std::stringstream ssF; ssF<<sim.path2file<<"/pressureIterStats.dat";
  std::ofstream pfile(ssF.str().c_str(), std::ofstream::app);
  if(sim.step==0) pfile<<"step time dt iter relDF"<<std::endl;
  pfile<<sim.step<<" "<<sim.time<<" "<<sim.dt<<" "<<iter<<" "<<relDF<<std::endl;
  }
  */

  check("IterativePressurePenalization");
}

CubismUP_3D_NAMESPACE_END
