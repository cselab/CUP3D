//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"

//#include "../poisson/PoissonSolverScalarACC_freespace.h"
//#include "../poisson/PoissonSolverScalarACC.h"

#include "../poisson/PoissonSolverUnbounded.h"
#include "../poisson/PoissonSolverPeriodic.h"
#include "../poisson/PoissonSolverMixed.h"

class OperatorGradP : public GenericLabOperator
{
 private:
  const double dt;
  const Real extent[3];

 public:
  OperatorGradP(double _dt,const Real ext[3]):dt(_dt),extent{ext[0],ext[1],ext[2]}
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2; stencil_end[1] = 2; stencil_end[2] = 2;
  }

  ~OperatorGradP() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // p contains the pressure correction after the Poisson solver
      const Real Uf = o(ix,iy,iz).u + fac*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p);
      const Real Vf = o(ix,iy,iz).v + fac*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p);
      const Real Wf = o(ix,iy,iz).w + fac*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p);
      const Real US=o(ix,iy,iz).tmpU, VS=o(ix,iy,iz).tmpV, WS=o(ix,iy,iz).tmpW;
      o(ix,iy,iz).u = ( Uf + o(ix,iy,iz).chi * US) / (1+o(ix,iy,iz).chi);
      o(ix,iy,iz).v = ( Vf + o(ix,iy,iz).chi * VS) / (1+o(ix,iy,iz).chi);
      o(ix,iy,iz).w = ( Wf + o(ix,iy,iz).chi * WS) / (1+o(ix,iy,iz).chi);
    }
  }
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
 protected:
  PoissonSolver * pressureSolver;

 public:
  CoordinatorPressure(SimulationData & s) : GenericCoordinator(s)
  {
    if(sim.bUseFourierBC)
    pressureSolver = new PoissonSolverPeriodic(sim);
    else if (sim.bUseUnboundedBC)
    pressureSolver = new PoissonSolverUnbounded(sim);
    else pressureSolver = new PoissonSolverMixed(sim);
    sim.pressureSolver = pressureSolver;
  }

  void operator()(const double dt)
  {

    pressureSolver->solve();

    { //pressure correction dudt* = - grad P / rho
      const int nthreads = omp_get_max_threads();
      std::vector<OperatorGradP*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i) diff[i] = new OperatorGradP(dt, sim.extent);

      compute<OperatorGradP>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    check("pressure - end");
  }

  std::string getName()
  {
    return "Pressure";
  }
};
#endif
