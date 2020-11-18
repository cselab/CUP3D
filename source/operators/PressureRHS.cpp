//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureRHS.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

//static inline PenalizationBlock* getPenalBlockPtr(
//  PenalizationGridMPI*const grid, const int blockID) {
//  assert(grid not_eq nullptr);
//  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();
//  return (PenalizationBlock*) vInfo[blockID].ptrBlock;
//}

template<bool implicitPenalization>
struct KernelPressureRHS : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  const Real dt = sim.dt;
  //PoissonSolver * const solver = sim.pressureSolver;
  PoissonSolverAMR * const solver = sim.pressureSolver;
  std::vector<int> & bElemTouchSurf;
  ObstacleVector * const obstacle_vector = sim.obstacle_vector;
  const int nShapes = obstacle_vector->nObstacles();
  // non-const non thread safe:
  std::vector<Real> sumRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> posRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> negRHS = std::vector<Real>(nShapes, 0);
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil =
  #ifdef PENAL_THEN_PRES
    nShapes>0 ? StencilInfo(-1,-1,-1, 2,2,2, false,
                                    {FE_U,FE_V,FE_W,FE_TMPU,FE_TMPV,FE_TMPW} ) :
  #endif
                StencilInfo(-1,-1,-1, 2,2,2, false, {FE_U, FE_V, FE_W});


  //FluxCorrectionMPI<FluxCorrection<BlockType,FluidElement>,FluidGridMPI> Corrector;

  KernelPressureRHS(SimulationData& s, std::vector<int> & bETS) :sim(s), bElemTouchSurf(bETS) 
  {
    //Corrector.prepare(*(sim.grid));
  }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {

    const Real h = info.h_gridpoint, fac = 0.5*h*h/dt;
    Real* __restrict__ const ret = solver->data.data() + solver->_offset(info);
    const unsigned SY=BlockType::sizeX, SZ=BlockType::sizeX*BlockType::sizeY;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      ret[SZ*iz +SY*iy +ix] = fac*(LE.u-LW.u + LN.v-LS.v + LB.w-LF.w);
    }

    //BlockCase<BlockType> * tempCase = (BlockCase<BlockType> *)(info.auxiliary);
    //BlockType::ElementType * faceXm = nullptr;
    //BlockType::ElementType * faceXp = nullptr;
    //BlockType::ElementType * faceYm = nullptr;
    //BlockType::ElementType * faceYp = nullptr;
    //BlockType::ElementType * faceZp = nullptr;
    //BlockType::ElementType * faceZm = nullptr;
    //if (tempCase != nullptr)
    //{
    //  faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
    //  faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
    //  faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
    //  faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
    //  faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
    //  faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
    //}
    //if (faceXm != nullptr)
    //{
    //  int ix = 0;
    //  for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    //  for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    //  {
    //    faceXm[iy + FluidBlock::sizeY * iz].clear();
    //    faceXm[iy + FluidBlock::sizeY * iz].p = fac *(lab(ix-1,iy,iz).u + lab(ix,iy,iz).u);
    //  }
    //}
    //if (faceXp != nullptr)
    //{
    //  int ix = FluidBlock::sizeX-1;
    //  for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    //  for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    //  {
    //    faceXp[iy + FluidBlock::sizeY * iz].clear();
    //    faceXp[iy + FluidBlock::sizeY * iz].p = - fac *(lab(ix+1,iy,iz).u + lab(ix,iy,iz).u);
    //  }
    //}   
    //if (faceYm != nullptr)
    //{
    //  int iy = 0;
    //  for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    //  for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    //  {
    //    faceYm[ix + FluidBlock::sizeX * iz].clear();
    //    faceYm[ix + FluidBlock::sizeX * iz].p = fac *(lab(ix,iy-1,iz).v + lab(ix,iy,iz).v);
    //  }
    //}
    //if (faceYp != nullptr)
    //{
    //  int iy = FluidBlock::sizeY-1;
    //  for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    //  for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    //  {
    //    faceYp[ix + FluidBlock::sizeX * iz].clear();
    //    faceYp[ix + FluidBlock::sizeX * iz].p = - fac *(lab(ix,iy+1,iz).v + lab(ix,iy,iz).v);
    //  }
    //}
    //if (faceZm != nullptr)
    //{
    //  int iz = 0;
    //  for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    //  for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    //  {
    //    faceZm[ix + FluidBlock::sizeX * iy].clear();
    //    faceZm[ix + FluidBlock::sizeX * iy].p = fac *(lab(ix,iy,iz-1).w + lab(ix,iy,iz).w);
    //  }
    //}
    //if (faceZp != nullptr)
    //{
    //  int iz = FluidBlock::sizeZ-1;
    //  for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    //  for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    //  {
    //    faceZp[ix + FluidBlock::sizeX * iy].clear();
    //    faceZp[ix + FluidBlock::sizeX * iy].p = - fac *(lab(ix,iy,iz+1).w + lab(ix,iy,iz).w);
    //  }
    //}

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
    if (obstblocks[info.blockID] == nullptr) return;

    const CHIMAT & __restrict__ CHI = obstblocks[info.blockID]->chi;
    const size_t offset = solver->_offset(info);
    Real* __restrict__ const ret = solver->data.data();
    const unsigned SY=BlockType::sizeX, SZ=BlockType::sizeX*BlockType::sizeY;
    const int obstID = obstacle->obstacleID; assert(obstID < nShapes);

    // Obstacle-specific lambda, useful for gradually adding an obstacle to the flow.
    const Real rampUp = obstacle->lambda_factor;
    // lambda = 1/dt hardcoded for expl time int, other options are wrong.
    const Real lamDt = rampUp * (implicitPenalization? sim.lambda * dt : 1.0);
    #ifdef PENAL_THEN_PRES
      const Real h = info.h_gridpoint, fac = 0.5*h*h/dt;
    #endif

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (lab(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      const size_t idx = offset + SZ * iz + SY * iy + ix;
      assert(idx < bElemTouchSurf.size());
      const Real X=CHI[iz][iy][ix], penalFac = not implicitPenalization? X*lamDt
                                               : X * lamDt/(1 + lamDt * X);
      #ifdef PENAL_THEN_PRES
        const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
        const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
        const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
        const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
        const Real srcBulk = - penalFac * fac * divUs;
      #else
        const Real srcBulk = - penalFac * ret[idx];
      #endif
      bElemTouchSurf[idx] = obstID;
      posRHS[obstID] += CHI[iz][iy][ix];
      sumRHS[obstID] += srcBulk;
      ret[idx] += srcBulk;
      negRHS[obstID] += std::fabs(ret[idx]);
    }
  }
};

struct KernelPressureRHS_nonUniform : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  const Real dt = sim.dt;
  //PoissonSolver * const solver = sim.pressureSolver;
  PoissonSolverAMR * const solver = sim.pressureSolver;
  std::vector<int> & bElemTouchSurf;
  ObstacleVector * const obstacle_vector = sim.obstacle_vector;
  const int nShapes = obstacle_vector->nObstacles();
  // non-const non thread safe:
  std::vector<Real> sumRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> posRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> negRHS = std::vector<Real>(nShapes, 0);
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_CHI,FE_U,FE_V,FE_W}};

  KernelPressureRHS_nonUniform(SimulationData & s, std::vector<int> & bETS) :
    sim(s), bElemTouchSurf(bETS) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    Real* __restrict__ const ret = solver->data.data() + solver->_offset(info);
    const unsigned SY=BlockType::sizeX, SZ=BlockType::sizeX*BlockType::sizeY;

    const BlkCoeffX &cx =o.fd_cx.first, &cy =o.fd_cy.first, &cz =o.fd_cz.first;
    const Real invdt = 1.0 / dt;

    for(unsigned iz=0; iz < (unsigned) FluidBlock::sizeZ; ++iz)
    for(unsigned iy=0; iy < (unsigned) FluidBlock::sizeY; ++iy)
    for(unsigned ix=0; ix < (unsigned) FluidBlock::sizeX; ++ix) {
      Real h[3]; info.spacing(h, ix, iy, iz);
      const FluidElement &L  = lab(ix  ,iy,  iz  );
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real dudx = __FD_2ND(ix, cx, LW.u, L.u, LE.u);
      const Real dvdy = __FD_2ND(iy, cy, LS.v, L.v, LN.v);
      const Real dwdz = __FD_2ND(iz, cz, LF.w, L.w, LB.w);
      ret[SZ*iz +SY*iy +ix] = h[0]*h[1]*h[2]*invdt * (dudx + dvdy + dwdz);
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
    if (obstblocks[info.blockID] == nullptr) return;

    const  CHIMAT & __restrict__  CHI = obstblocks[info.blockID]->chi;
    const FluidBlock & b = * (FluidBlock *) info.ptrBlock;

    const size_t offset = solver->_offset(info);
    Real* __restrict__ const ret = solver->data.data();
    const unsigned SY=BlockType::sizeX, SZ=BlockType::sizeX*BlockType::sizeY;

    const BlkCoeffX &cx =b.fd_cx.first, &cy =b.fd_cy.first, &cz =b.fd_cz.first;
    const int obstID = obstacle->obstacleID; assert(obstID < nShapes);

    for(unsigned iz=0; iz < (unsigned) FluidBlock::sizeZ; ++iz)
    for(unsigned iy=0; iy < (unsigned) FluidBlock::sizeY; ++iy)
    for(unsigned ix=0; ix < (unsigned) FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L = lab(ix,iy,iz);
      if (L.chi > CHI[iz][iy][ix]) continue;
      const size_t idx = offset + SZ * iz + SY * iy + ix;
      assert(idx < bElemTouchSurf.size());

      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real dXdx = __FD_2ND(ix, cx, LW.chi, L.chi, LE.chi);
      const Real dXdy = __FD_2ND(iy, cy, LS.chi, L.chi, LN.chi);
      const Real dXdz = __FD_2ND(iz, cz, LF.chi, L.chi, LB.chi);

      const bool bSurface = dXdx * dXdx + dXdy * dXdy + dXdz * dXdz > 0;
      if(bSurface) bElemTouchSurf[idx] = obstID;
      const Real srcBulk = - CHI[iz][iy][ix] * ret[idx];
      sumRHS[obstID] += srcBulk;
      ret[idx] += srcBulk;
      posRHS[obstID] += bSurface * ret[idx] * (ret[idx] > 0);
      negRHS[obstID] -= bSurface * ret[idx] * (ret[idx] < 0);
    }
  }
};

struct KernelFinalizePerimeters : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();
  //PoissonSolver * const solver;
  PoissonSolverAMR * const solver;
  const std::vector<int> & bElemTouchSurf;
  const std::vector<Real> & corrFactors;

  KernelFinalizePerimeters(FluidGridMPI* const g, PoissonSolverAMR* ps, /* PoissonSolver* ps, */
    std::vector<int> & bETS, const std::vector<Real> & corr): grid(g),
    solver(ps), bElemTouchSurf(bETS), corrFactors(corr) {}

  void visit(Obstacle* const obstacle)
  {
    const unsigned SY=BlockType::sizeX, SZ=BlockType::sizeX*BlockType::sizeY;
    const int obstID = obstacle->obstacleID;
    const Real corr = corrFactors[obstID];

    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        const cubism::BlockInfo& info = vInfo[i];
        if(obstblocks[info.blockID] == nullptr) continue;

        //const CHIMAT & __restrict__ CHI = obstblocks[info.blockID]->chi;
        const size_t offset = solver->_offset(info);
        Real* __restrict__ const ret = solver->data.data();

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
          const size_t idx = offset + SZ * iz + SY * iy + ix;
          if (bElemTouchSurf[idx] not_eq obstID) continue;
          //ret[idx] -= corr * CHI[iz][iy][ix];
          ret[idx] -= corr * std::fabs(ret[idx]);
        }
      }
    }
  }
};

struct PressureRHSObstacleVisitor : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
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
        const CHIMAT & __restrict__ CHI = pos->chi;
        //const CHIMAT & __restrict__ SDF = pos->sdf;

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

}

PressureRHS::PressureRHS(SimulationData & s) : Operator(s)
{
  // no need to allocate stretched mesh stuff here because we will never read
  // grid spacing from this grid!
  if(sim.rank==0) printf("Allocating the penalization helper grid.\n");

  std::cout << "penalizationGrid not allocated." << std::endl;
  //penalizationGrid = new PenalizationGridMPI(
  //  1, //these arguments are not used by Cubism-AMR
  //  1, //these arguments are not used by Cubism-AMR
  //  1, //these arguments are not used by Cubism-AMR
  //  sim.bpdx,
  //  sim.bpdy,
  //  sim.bpdz, sim.maxextent,sim.levelStart,sim.levelMax,sim.app_comm,
  //            (sim.BCx_flag == periodic),
  //            (sim.BCy_flag == periodic),
  //            (sim.BCz_flag == periodic));
}

//PressureRHS::~PressureRHS() { /*delete penalizationGrid;*/ }


#define PRESRHS_LOOP(T) do {                                                 \
      std::vector< T *> K(nthreads, nullptr);                                \
      for(int i=0;i<nthreads;++i) K[i] = new T (sim,elemTouchSurf);          \
      compute< T >(K);                                                       \
      for(size_t j = 0; j<nShapes; ++j) {                                    \
        for(int i=0; i<nthreads; ++i) sumRHS[j] += K[i]->sumRHS[j];          \
        for(int i=0; i<nthreads; ++i) posRHS[j] += K[i]->posRHS[j];          \
        for(int i=0; i<nthreads; ++i) negRHS[j] += K[i]->negRHS[j];          \
      }                                                                      \
      for(int i=0; i<nthreads; i++) delete K[i];                             \
    } while(0)

void PressureRHS::operator()(const double dt)
{
  //place onto p: ( div u^(t+1) - div u^* ) / dt
  //where i want div u^(t+1) to be equal to div udef
  sim.pressureSolver->reset();

  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();
  const size_t nShapes = sim.obstacle_vector->nObstacles();
  std::vector<Real> corrFactors(nShapes, 0);
  const int nthreads = omp_get_max_threads();
  std::vector<double> sumRHS(nShapes,0), posRHS(nShapes,0), negRHS(nShapes,0);
  std::vector<int> elemTouchSurf(sim.pressureSolver->data.size(), -1);

  #ifdef PENAL_THEN_PRES
    sim.startProfiler("PresRHS Udef");
    if(nShapes > 0)
    { //zero fields, going to contain Udef:
      #pragma omp parallel for schedule(static)
      for(size_t i=0; i<vInfo.size(); i++) {
        FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
          b(ix,iy,iz).tmpU = 0; b(ix,iy,iz).tmpV = 0; b(ix,iy,iz).tmpW = 0;
        }
      }
      //store deformation velocities onto tmp fields:
      ObstacleVisitor* visitor = new PressureRHSObstacleVisitor(grid);
      sim.obstacle_vector->Accept(visitor);
      delete visitor;
    }
    sim.stopProfiler();
  #endif

  sim.startProfiler("PresRHS Kernel");
  if(sim.bUseStretchedGrid) PRESRHS_LOOP(KernelPressureRHS_nonUniform);
  else {
    if(sim.bImplicitPenalization) PRESRHS_LOOP(KernelPressureRHS<1>);
    else PRESRHS_LOOP(KernelPressureRHS<0>);
  }
  sim.stopProfiler();

  if(nShapes == 0) return; // no need to deal with obstacles perimeters
  // non-divergence free obstacles may have int div u not_eq 0
  // meaning that they will have net out/in flow
  // usually it is a small number and here we correct this:

  const auto& COMM = sim.app_comm;
  MPI_Allreduce(MPI_IN_PLACE, sumRHS.data(), nShapes, MPI_DOUBLE,MPI_SUM, COMM);
  MPI_Allreduce(MPI_IN_PLACE, posRHS.data(), nShapes, MPI_DOUBLE,MPI_SUM, COMM);
  MPI_Allreduce(MPI_IN_PLACE, negRHS.data(), nShapes, MPI_DOUBLE,MPI_SUM, COMM);
  for(size_t j = 0; j<nShapes; ++j) {
      corrFactors[j] = sumRHS[j] / std::max(negRHS[j], (double) EPS);
      //corrFactors[j] = sumRHS[j] / std::max(posRHS[j], (double) EPS);
  }

  {
    sim.startProfiler("PresRHS Correct");
    ObstacleVisitor*K = new KernelFinalizePerimeters(
      grid, sim.pressureSolver, elemTouchSurf, corrFactors);
    sim.obstacle_vector->Accept(K); // accept you son of a french cow
    delete K;
    sim.stopProfiler();
  }

  check("PressureRHS");
}

CubismUP_3D_NAMESPACE_END
