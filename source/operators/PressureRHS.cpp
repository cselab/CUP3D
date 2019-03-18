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

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

template<int withObstacles = 1>
class KernelPressureRHS
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
    return withObstacles ? divUf - L.chi*divUs : divUf;
  }

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,6,1,2,3,5,6,7);


  KernelPressureRHS(double _dt, double lambda, const Real B[3], const std::array<Real, 3> &E,
    PoissonSolver* ps) : dt(_dt), lamdt(_dt*lambda), fadeLen{B[0],B[1],B[2]},
    ext{E[0],E[1],E[2]}, iFade{1/(B[0]+EPS),1/(B[1]+EPS),1/(B[2]+EPS)}, solver(ps) {}
  ~KernelPressureRHS() {}

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

template<int withObstacles = 1>
class KernelPressureRHS_nonUniform
{
 private:
  Real dt, invdt, lamdt;
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

  inline Real RHS(Lab&l, const int ix, const int iy, const int iz,
    const BlkCoeffX &cx, const BlkCoeffX &cy, const BlkCoeffX &cz) const
  {
    const FluidElement& L  = l(ix,  iy,  iz);
    const FluidElement& LW = l(ix-1,iy,  iz  ), & LE = l(ix+1,iy,  iz  );
    const FluidElement& LS = l(ix,  iy-1,iz  ), & LN = l(ix,  iy+1,iz  );
    const FluidElement& LF = l(ix,  iy,  iz-1), & LB = l(ix,  iy,  iz+1);
    const Real dudx = __FD_2ND(ix, cx, LW.u, L.u, LE.u);
    const Real dvdy = __FD_2ND(iy, cy, LS.v, L.v, LN.v);
    const Real dwdz = __FD_2ND(iz, cz, LF.w, L.w, LB.w);
    const Real divUf = dudx + dvdy + dwdz;
    const Real dusdx = __FD_2ND(ix, cx, LW.tmpU, L.tmpU, LE.tmpU);
    const Real dvsdy = __FD_2ND(iy, cy, LS.tmpV, L.tmpV, LN.tmpV);
    const Real dwsdz = __FD_2ND(iz, cz, LF.tmpW, L.tmpW, LB.tmpW);
    const Real divUs = dusdx + dvsdy + dwsdz;
    return withObstacles ? divUf - L.chi*divUs : divUf;
  }

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1,2,2,2,false, 6, 1,2,3,5,6,7);


  KernelPressureRHS_nonUniform(double _dt, double lambda, const Real buf[3],
  const std::array<Real, 3> &E, PoissonSolver* ps): dt(_dt), invdt(1/_dt), lamdt(_dt*lambda),
  fadeLen{buf[0],buf[1],buf[2]}, ext{E[0],E[1],E[2]}, iFade{1/(buf[0]+EPS), 1/(buf[1]+EPS), 1/(buf[2]+EPS)}, solver(ps) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const size_t offset = solver->_offset_ext(info);
    // FD coefficients for first derivative
    const BlkCoeffX &cx =o.fd_cx.first, &cy =o.fd_cy.first, &cz =o.fd_cz.first;
    if( not _is_touching(o) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        Real h[3]; info.spacing(h, ix, iy, iz);
        const Real fac = h[0]*h[1]*h[2]*invdt;
        const Real RHS_ = RHS(lab, ix,iy,iz, cx,cy,cz);
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS_ );
        //o(ix,iy,iz).p = fac * RHS(lab, ix,iy,iz, pFac); //will break t>0
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        Real h[3]; info.spacing(h, ix, iy, iz);
        const Real fac = h[0]*h[1]*h[2]*invdt;
        const Real RHS_ = fade(info, ix,iy,iz) * RHS(lab, ix,iy,iz, cx,cy,cz);
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS_ );
        //o(ix,iy,iz).p = fac * RHS_; //will break t>0
      }
    }
  }
};

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

void PressureRHS::operator()(const double dt)
{
  sim.startProfiler("PresRHS Udef");
  if(sim.obstacle_vector->nObstacles() > 0)
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

  sim.startProfiler("PresRHS Kernel");
  //place onto p: ( div u^(t+1) - div u^* ) / dt
  //where i want div u^(t+1) to be equal to div udef
  if(sim.bUseStretchedGrid)
  {
    if(sim.obstacle_vector->nObstacles())
    {
      const KernelPressureRHS_nonUniform<1> K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
      compute<KernelPressureRHS_nonUniform<1>>(K);
    }
    else
    {
      const KernelPressureRHS_nonUniform<0> K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
      compute<KernelPressureRHS_nonUniform<0>>(K);
    }
  }
  else
  {
    if(sim.obstacle_vector->nObstacles())
    {
      const KernelPressureRHS<1> K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
      compute<KernelPressureRHS<1>>(K);
    }
    else
    {
      const KernelPressureRHS<0> K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
      compute<KernelPressureRHS<0>>(K);
    }
  }


  sim.stopProfiler();

  check("PressureRHS");
}

CubismUP_3D_NAMESPACE_END
