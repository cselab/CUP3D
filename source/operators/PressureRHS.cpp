//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/PressureRHS.h"
#include "operators/PenalizationObstacleVisitor.h"
#include "poisson/PoissonSolver.h"

CubismUP_3D_NAMESPACE_BEGIN

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

  inline Real RHS(Lab&l, const int x,const int y,const int z,const Real F) const
  {
    const FluidElement & L  = l(x,  y,  z);
    const FluidElement & LW = l(x-1,y,  z  ), & LE = l(x+1,y,  z  );
    const FluidElement & LS = l(x,  y-1,z  ), & LN = l(x,  y+1,z  );
    const FluidElement & LF = l(x,  y,  z-1), & LB = l(x,  y,  z+1);
    const Real divUt = LE.u-LW.u + LN.v-LS.v + LB.w-LF.w;
    #if PENAL_TYPE==2
      const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
      const Real dXx_dux = (LE.chi-LW.chi)*( L.u - F * (LE.p-LW.p) - L.tmpU );
      const Real dXy_duy = (LN.chi-LS.chi)*( L.v - F * (LN.p-LS.p) - L.tmpV );
      const Real dXz_duz = (LB.chi-LF.chi)*( L.w - F * (LB.p-LF.p) - L.tmpW );
      const Real X = L.chi, fac1= lamdt*(X -X*X) -X, fac2= -lamdt/(1+lamdt*X);
      return divUt + fac1*divUs + fac2*(dXx_dux + dXy_duy + dXz_duz);
    #elif PENAL_TYPE==1
      const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
      return divUt - L.chi*divUs;
    #else
      const Real uFx_dXx = ( LE.chi - LW.chi ) * L.tmpU;
      const Real uFy_dXy = ( LN.chi - LS.chi ) * L.tmpV;
      const Real uFz_dXz = ( LB.chi - LF.chi ) * L.tmpW;
      return divUt + uFx_dXx + uFy_dXy + uFz_dXz;
    #endif
  }

 public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  #if PENAL_TYPE==2
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,8,0,1,2,3,4,5,6,7);
  #elif PENAL_TYPE==1
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,6,1,2,3,5,6,7);
  #else
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,4,0,1,2,3);
  #endif


  KernelPressureRHS(double _dt, double lambda, const Real B[3], const Real E[3],
    PoissonSolver* ps) : dt(_dt), lamdt(_dt*lambda), fadeLen{B[0],B[1],B[2]},
    ext{E[0],E[1],E[2]}, iFade{1/(B[0]+EPS),1/(B[1]+EPS),1/(B[2]+EPS)}, solver(ps) {}
  ~KernelPressureRHS() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt, pFac = .5*dt/h;
    const size_t offset = solver->_offset_ext(info);
    if( not _is_touching(o) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS(lab, ix,iy,iz, pFac) );
        //o(ix,iy,iz).p = fac * RHS(lab, ix,iy,iz, pFac); //will break t>0
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        const Real RHS_ = fade(info, ix,iy,iz) * RHS(lab, ix,iy,iz, pFac);
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS_);
        //o(ix,iy,iz).p = fac * RHS_; //will break t>0
      }
    }
  }
};

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
    const Real divUt = dudx + dvdy + dwdz;
    #if PENAL_TYPE!=0
      const Real dusdx = __FD_2ND(ix, cx, LW.tmpU, L.tmpU, LE.tmpU);
      const Real dvsdy = __FD_2ND(iy, cy, LS.tmpV, L.tmpV, LN.tmpV);
      const Real dwsdz = __FD_2ND(iz, cz, LF.tmpW, L.tmpW, LB.tmpW);
      const Real divUs = dusdx + dvsdy + dwsdz;
    #endif
    #if PENAL_TYPE!=1
      const Real dXdx = __FD_2ND(ix, cx, LW.chi, L.chi, LE.chi);
      const Real dXdy = __FD_2ND(iy, cy, LS.chi, L.chi, LN.chi);
      const Real dXdz = __FD_2ND(iz, cz, LF.chi, L.chi, LB.chi);
    #endif

    #if PENAL_TYPE==2
      const Real dpdx = __FD_2ND(ix, cx, LW.p, L.p, LE.p);
      const Real dpdy = __FD_2ND(iy, cy, LS.p, L.p, LN.p);
      const Real dpdz = __FD_2ND(iz, cz, LF.p, L.p, LB.p);
      const Real dXx_dux = dXdx * ( L.u - dt * dpdx - L.tmpU );
      const Real dXy_duy = dXdy * ( L.v - dt * dpdy - L.tmpV );
      const Real dXz_duz = dXdz * ( L.w - dt * dpdz - L.tmpW );
      //this assumes div(u^t+1)= chi div(u_s), underestimates sphere Cd:
      const Real X = L.chi, fac1= lamdt*(X -X*X) -X, fac2= -lamdt/(1+lamdt*X);
      return divUt + fac1*divUs + fac2*(dXx_dux + dXy_duy + dXz_duz);
    #elif PENAL_TYPE==1
      return divUt - L.chi*divUs;
    #else
      const Real uFx_dXx=dXdx*L.tmpU, uFy_dXy=dXdy*L.tmpV, uFz_dXz=dXdz*L.tmpW;
      return divUt + uFx_dXx + uFy_dXy + uFz_dXz;
    #endif
  }

 public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  #if PENAL_TYPE==2
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,8,0,1,2,3,4,5,6,7);
  #elif PENAL_TYPE==1
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,6,1,2,3,5,6,7);
  #else
  const StencilInfo stencil=StencilInfo(-1,-1,-1,2,2,2,false,4,0,1,2,3);
  #endif


  KernelPressureRHS_nonUniform(double _dt, double lambda, const Real buf[3],
  const Real E[3], PoissonSolver* ps): dt(_dt), invdt(1/_dt), lamdt(_dt*lambda),
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
        solver->_cub2fftw(offset, iz,iy,ix, fac * fac );
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
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS_);
        //o(ix,iy,iz).p = fac * RHS_; //will break t>0
      }
    }
  }
};

void PressureRHS::operator()(const double dt)
{
  #if PENAL_TYPE!=0
  sim.startProfiler("PresRHS Uobst.");
  { //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<vInfo.size(); i++) {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
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
  sim.stopProfiler();
  #endif

  sim.startProfiler("PresRHS Kernel");
  //place onto p: ( div u^(t+1) - div u^* ) / dt
  //where i want div u^(t+1) to be equal to div udef
  const KernelPressureRHS K(dt, sim.lambda, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);
  compute<KernelPressureRHS>(K);
  sim.stopProfiler();

  check("pressure rhs - end");
}

CubismUP_3D_NAMESPACE_END
