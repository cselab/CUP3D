//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/PressureRHS.h"
#include "poisson/PoissonSolver.h"
#include "obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

static inline PenalizationBlock* getPenalBlockPtr(
  PenalizationGridMPI*const grid, const int blockID) {
  assert(grid not_eq nullptr);
  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();
  return (PenalizationBlock*) vInfo[blockID].ptrBlock;
}

struct KernelPressureRHS : public ObstacleVisitor
{
  const Real dt;
  PoissonSolver * const solver;
  PenalizationGridMPI * const penGrid;
  ObstacleVector * const obstacle_vector;
  const size_t nShapes = obstacle_vector->nObstacles();
  // non-const non thread safe:
  Real * const sumRHS = (Real*) calloc(nShapes, sizeof(Real));
  Real * const absRHS = (Real*) calloc(nShapes, sizeof(Real));
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 4, 0,1,2,3);

  KernelPressureRHS(PoissonSolver* pois, PenalizationGridMPI* penG,
    ObstacleVector* oVec, double _dt) : dt(_dt), solver(pois), penGrid(penG),
    obstacle_vector(oVec) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt;
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      ret[SZ*iz +SY*iy +SX*ix] = fac*(LE.u-LW.u + LN.v-LS.v + LB.w-LF.w);
    }

    if(nShapes == 0) return; // no need to account for obstacles

    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr && lab_ptr == nullptr);
    info_ptr =  & info; lab_ptr =   & lab;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr; lab_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    assert(info_ptr not_eq nullptr && lab_ptr not_eq nullptr);
    const BlockInfo& info = * info_ptr;
    Lab& lab = * lab_ptr;
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const  CHIMAT & __restrict__  CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    const Real h = info.h_gridpoint, fac = 0.5 * h * h / dt;
    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();

    PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;
    const size_t obstID = obstacle->obstacleID; assert(obstID < nShapes);

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L = lab(ix,iy,iz); t(ix,iy,iz).uPres = L.chi;
      if (CHI[iz][iy][ix] <= 0 || L.chi > CHI[iz][iy][ix]) continue;
      Real p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];
      const Real U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };

      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real dUdChiX = ( LE.chi - LW.chi ) * ( U_TOT[0] - L.u );
      const Real dVdChiY = ( LN.chi - LS.chi ) * ( U_TOT[1] - L.v );
      const Real dWdChiZ = ( LB.chi - LF.chi ) * ( U_TOT[2] - L.w );
      const Real srcBulk = - fac * L.chi * (LE.u-LW.u + LN.v-LS.v + LB.w-LF.w);
      const Real srcPerimeter = fac * ( dUdChiX + dVdChiY + dWdChiZ );

      sumRHS[obstID]  += srcPerimeter + srcBulk;
      t(ix,iy,iz).rhs0 = std::fabs( srcPerimeter );
      absRHS[obstID]  += std::fabs( srcPerimeter );
      ret[SZ*iz + SY*iy + SX*ix] += srcPerimeter + srcBulk;
    }
  }
};

struct KernelPressureRHS_nonUniform : public ObstacleVisitor
{
  const Real dt, invdt = 1.0/dt;
  PoissonSolver * const solver;
  PenalizationGridMPI * const penGrid;
  ObstacleVector * const obstacle_vector;
  const size_t nShapes = obstacle_vector->nObstacles();
  // non-const non thread safe:
  Real * const sumRHS = (Real*) calloc(nShapes, sizeof(Real));
  Real * const absRHS = (Real*) calloc(nShapes, sizeof(Real));
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 4, 0,1,2,3);

  KernelPressureRHS_nonUniform(PoissonSolver* pois, PenalizationGridMPI* penG,
    ObstacleVector* oVec, double _dt) : dt(_dt), solver(pois), penGrid(penG),
    obstacle_vector(oVec) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;
    const BlkCoeffX &cx =o.fd_cx.first, &cy =o.fd_cy.first, &cz =o.fd_cz.first;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      Real h[3]; info.spacing(h, ix, iy, iz);
      const FluidElement &L  = lab(ix  ,iy,  iz  );
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real dudx = __FD_2ND(ix, cx, LW.u, L.u, LE.u);
      const Real dvdy = __FD_2ND(iy, cy, LS.v, L.v, LN.v);
      const Real dwdz = __FD_2ND(iz, cz, LF.w, L.w, LB.w);
      ret[SZ*iz +SY*iy +SX*ix] = h[0]*h[1]*h[2]*invdt * (dudx + dvdy + dwdz);
    }

    if(nShapes == 0) return; // no need to account for obstacles

    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr && lab_ptr == nullptr);
    info_ptr = & info;
    lab_ptr = & lab;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
    lab_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    assert(info_ptr not_eq nullptr && lab_ptr not_eq nullptr);
    const BlockInfo& info = * info_ptr;
    Lab& lab = * lab_ptr;
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const  CHIMAT & __restrict__  CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();

    const FluidBlock & b = * (FluidBlock *) info.ptrBlock;
    PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;
    const BlkCoeffX &cx =b.fd_cx.first, &cy =b.fd_cy.first, &cz =b.fd_cz.first;
    const size_t obstID = obstacle->obstacleID; assert(obstID < nShapes);

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L = lab(ix,iy,iz); t(ix,iy,iz).uPres = L.chi;
      if (CHI[iz][iy][ix] <= 0 || L.chi > CHI[iz][iy][ix]) continue;
      Real p[3]; info.pos(p, ix, iy, iz);
      Real h[3]; info.spacing(h, ix, iy, iz);
      const Real fac = h[0] * h[1] * h[2] * invdt;
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];
      const Real U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };

      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real dchidx = __FD_2ND(ix, cx, LW.chi, L.chi, LE.chi);
      const Real dchidy = __FD_2ND(iy, cy, LS.chi, L.chi, LN.chi);
      const Real dchidz = __FD_2ND(iz, cz, LF.chi, L.chi, LB.chi);
      const Real dudx = __FD_2ND(ix, cx, LW.u, L.u, LE.u);
      const Real dvdy = __FD_2ND(iy, cy, LS.v, L.v, LN.v);
      const Real dwdz = __FD_2ND(iz, cz, LF.w, L.w, LB.w);
      const Real dUdChiX = dchidx * ( U_TOT[0] - L.u );
      const Real dVdChiY = dchidy * ( U_TOT[1] - L.v );
      const Real dWdChiZ = dchidz * ( U_TOT[2] - L.w );
      const Real srcBulk = - fac * L.chi * (dudx + dvdy + dwdz);
      const Real srcPerimeter = fac * ( dUdChiX + dVdChiY + dWdChiZ );

      sumRHS[obstID]  += srcPerimeter + srcBulk;
      t(ix,iy,iz).rhs0 = std::fabs( srcPerimeter );
      absRHS[obstID]  += std::fabs( srcPerimeter );
      ret[SZ*iz + SY*iy + SX*ix] += srcPerimeter + srcBulk;
    }
  }
};

struct KernelFinalizePerimeters : public ObstacleVisitor
{
  PoissonSolver * const solver;
  PenalizationGridMPI * const penGrid;
  const std::vector<cubism::BlockInfo>& vInfo = penGrid->getBlocksInfo();
  const std::vector<Real> & corrFactors;

  KernelFinalizePerimeters(PoissonSolver* ps, PenalizationGridMPI* pen,
   const std::vector<Real>&corr): solver(ps), penGrid(pen), corrFactors(corr) {}

  void visit(Obstacle* const obstacle)
  {
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;
    const size_t obstID = obstacle->obstacleID;

    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        const cubism::BlockInfo& info = vInfo[i];
        const ObstacleBlock*const o = obstblocks[info.blockID];
        if(o == nullptr) continue;

        const CHIMAT & __restrict__ CHI = o->chi;
        const PenalizationBlock& t = * getPenalBlockPtr(penGrid, info.blockID);
        Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
          // wherever relevant, chi was stored in t(ix,iy,iz).uPres
          // and rhs values at object boundaries in t(ix,iy,iz).rhs0
          if (CHI[iz][iy][ix]<=0 || t(ix,iy,iz).uPres>CHI[iz][iy][ix]) continue;
          ret[SZ*iz + SY*iy + SX*ix] -= corrFactors[obstID] * t(ix,iy,iz).rhs0;
        }
      }
    }
  }
};

}

PressureRHS::PressureRHS(SimulationData & s) : Operator(s)
{
  // no need to allocate stretched mesh stuff here because we will never read
  // grid spacing from this grid!
  if(sim.rank==0) printf("Allocating the penalization helper grid.\n");

  penalizationGrid = new PenalizationGridMPI(
    sim.nprocsx, sim.nprocsy, sim.nprocsz, sim.local_bpdx,
    sim.local_bpdy, sim.local_bpdz, sim.maxextent, sim.app_comm);
}

PressureRHS::~PressureRHS() { delete penalizationGrid; }

void PressureRHS::operator()(const double dt)
{
  //place onto p: ( div u^(t+1) - div u^* ) / dt
  //where i want div u^(t+1) to be equal to div udef
  sim.pressureSolver->reset();

  const size_t nShapes = sim.obstacle_vector->nObstacles();
  std::vector<Real> corrFactors(nShapes, 0);
  const int nthreads = omp_get_max_threads();
  std::vector<double> sumRHS(nShapes, 0), absRHS(nShapes, 0);

  sim.startProfiler("PresRHS Kernel");
  if(sim.bUseStretchedGrid)
  {
    std::vector<KernelPressureRHS_nonUniform*> K(nthreads, nullptr);
    #pragma omp parallel for schedule(static, 1)
    for(int i=0; i<nthreads; ++i)
      K[i] = new KernelPressureRHS_nonUniform(sim.pressureSolver,
        penalizationGrid, sim.obstacle_vector, dt);

    compute<KernelPressureRHS_nonUniform>(K);

    for(size_t j = 0; j<nShapes; ++j) {
      for(int i=0; i<nthreads; ++i) absRHS[j] += K[i]->absRHS[j];
      for(int i=0; i<nthreads; ++i) sumRHS[j] += K[i]->sumRHS[j];
    }
    for(int i=0; i<nthreads; i++) delete K[i];
  }
  else
  {
    std::vector<KernelPressureRHS*> K(nthreads, nullptr);
    #pragma omp parallel for schedule(static, 1)
    for(int i=0; i<nthreads; ++i)
      K[i] = new KernelPressureRHS(sim.pressureSolver, penalizationGrid,
        sim.obstacle_vector, dt);

    compute<KernelPressureRHS>(K);

    for(size_t j = 0; j<nShapes; ++j) {
      for(int i=0; i<nthreads; ++i) absRHS[j] += K[i]->absRHS[j];
      for(int i=0; i<nthreads; ++i) sumRHS[j] += K[i]->sumRHS[j];
    }
    for(int i=0; i<nthreads; i++) delete K[i];
  }
  sim.stopProfiler();

  if(nShapes == 0) return; // no need to deal with obstacles perimeters

  const auto& COMM = sim.app_comm;
  MPI_Allreduce(MPI_IN_PLACE, sumRHS.data(), nShapes, MPI_DOUBLE,MPI_SUM, COMM);
  MPI_Allreduce(MPI_IN_PLACE, absRHS.data(), nShapes, MPI_DOUBLE,MPI_SUM, COMM);
  for(size_t j = 0; j<nShapes; ++j)
    corrFactors[j] = sumRHS[j] / std::max(absRHS[j], EPS);

  //if(0) // TODO DEBUG
  {
    sim.startProfiler("PresRHS Correct");
    ObstacleVisitor*K = new KernelFinalizePerimeters(sim.pressureSolver, penalizationGrid, corrFactors);
    sim.obstacle_vector->Accept(K); // accept you son of a french cow
    delete K;
    sim.stopProfiler();
  }

  check("PressureRHS");
}

CubismUP_3D_NAMESPACE_END
