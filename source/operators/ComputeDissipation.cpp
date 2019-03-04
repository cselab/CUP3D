//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Sid Verma in May 2018.
//

#include "operators/ComputeDissipation.h"
#include "utils/BufferedLogger.h"

CubismUP_3D_NAMESPACE_BEGIN

class KernelDissipation
{
  public:
  Real viscous = 0, press = 0, kinetic = 0;
  private:
  const double dt, nu;
  const Real extent[3];
  public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 4, 1,2,3,4);

  KernelDissipation(double _dt, const Real ext[3], Real _nu)
  : dt(_dt), nu(_nu), extent{ext[0],ext[1],ext[2]} { }
  ~KernelDissipation() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real h = info.h_gridpoint;
    const Real hCube = std::pow(h,3);
    const Real factor = 0.5/h/dt, invHsqr = 1/(h*h), dissipFactor = 2*hCube*nu;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);

      //shear stresses
      const Real D11 = factor*(LE.u - LW.u);
      const Real D22 = factor*(LN.v - LS.v);
      const Real D33 = factor*(LB.w - LF.w);
      const Real D12 = factor*(LN.u - LS.u + LE.v - LW.v)/2;
      const Real D13 = factor*(LB.u - LF.u + LE.w - LW.w)/2;
      const Real D23 = factor*(LN.w - LS.w + LB.v - LF.v)/2;

      // need to multiply this by 2*nu*h^3:
      const Real Sij = (D11*D11+D22*D22+D33*D33 + 2*(D12*D12+D13*D13+D23*D23));
      const Real dPdx = factor*(LE.p - LW.p);
      const Real dPdy = factor*(LN.p - LS.p);
      const Real dPdz = factor*(LB.p - LF.p);

      // Need to multiply this by h^3:
      const Real pressTerm = -(dPdx*L.u + dPdy*L.v + dPdz*L.w);

      const Real lapU_1 = invHsqr*(LE.u +LW.u +LN.u +LS.u +LB.u +LF.u - 6*L.u);
      const Real lapU_2 = invHsqr*(LE.v +LW.v +LN.v +LS.v +LB.v +LF.v - 6*L.v);
      const Real lapU_3 = invHsqr*(LE.w +LW.w +LN.w +LS.w +LB.w +LF.w - 6*L.w);
      // Store nu * u \cdot (nabla^2 u)
      // WARNING: multiply here by 0.5 coz later will multiply by dissipFactor=2*nu*hCube

      // need to multiply this by 2*nu*h^3:
      const Real laplacianTerm= (L.u*lapU_1 + L.v*lapU_2 + L.w*lapU_3)/2;
      viscous += (1-L.chi) * dissipFactor * (laplacianTerm + Sij);
      press += (1-L.chi) * hCube * pressTerm;
      kinetic += (1-L.chi) * hCube * (L.u*L.u + L.v*L.v + L.w*L.w)/2;
    }
  }
};

void ComputeDissipation::operator()(const double dt)
{
  sim.startProfiler("Dissip Kernel");
  const int nthreads = omp_get_max_threads();
  std::vector<KernelDissipation*> diss(nthreads, nullptr);
  #pragma omp parallel for schedule(static, 1)
  for(int i=0; i<nthreads; ++i)
    diss[i] = new KernelDissipation(dt, sim.extent, sim.nu);

  compute<KernelDissipation>(diss);
  sim.stopProfiler();

  sim.startProfiler("Dissip Reduce");
  double viscous=0.0, press=0.0, kinetic=0.0;
  for(int i=0; i<nthreads; i++)
  {
    kinetic += diss[i]->kinetic;
    viscous += diss[i]->viscous;
    press += diss[i]->press;
    delete diss[i];
  }

  double localSum[3] = {viscous, press, kinetic};
  double globalSum[3] = {0.0, 0.0, 0.0};
  MPI_Allreduce(localSum,globalSum,3, MPI_DOUBLE,MPI_SUM,grid->getCartComm());

  if(sim.rank==0)
  {
    std::stringstream &fileDissip = logger.get_stream("wakeDissipation.dat");
    if(sim.step==0)
    fileDissip<<"step_id time viscousTerm pressureTerm kineticEn"<<std::endl;

    fileDissip<<sim.step<<" "<<sim.time<<" "<<globalSum[0]<<" "
              <<globalSum[1] <<" "<<globalSum[2]<<std::endl;
  }
  sim.stopProfiler();

  check("dissipation - end");
}

CubismUP_3D_NAMESPACE_END
