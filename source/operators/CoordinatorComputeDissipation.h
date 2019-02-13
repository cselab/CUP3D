//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Sid Verma in May 2018.
//

#ifndef CubismUP_3D_CoordinatorDissipation_h
#define CubismUP_3D_CoordinatorDissipation_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "../utils/BufferedLogger.h"

class OperatorDissipation : public GenericLabOperator
{
  private:
  const double dt, nu;
  const Real extent[3];

  public:
  OperatorDissipation(double _dt, const Real ext[3], Real _nu)
  : dt(_dt), nu(_nu), extent{ext[0],ext[1],ext[2]}
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 4, 1,2,3,4);
    stencil_start[0] = -1;
    stencil_start[1] = -1;
    stencil_start[2] = -1;
    stencil_end[0] = 2;
    stencil_end[1] = 2;
    stencil_end[2] = 2;
  }
  ~OperatorDissipation() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint;
    //const Real hCube = pow(h,3);
    const Real factor = 0.5/h/dt;
    const Real invHsqr = 1.0/(h*h);
    //const Real dissipFactor = 2.0*hCube*nu;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {

      const Real chi   = lab(ix,iy,iz).chi;
      const Real uGrid = lab(ix,iy,iz).u;
      const Real vGrid = lab(ix,iy,iz).v;
      const Real wGrid = lab(ix,iy,iz).w;

      //shear stresses
      const Real D11 =    factor*(lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u);
      const Real D22 =    factor*(lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v);
      const Real D33 =    factor*(lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w);

      const Real D12 = .5*factor*(lab(ix,iy+1,iz).u - lab(ix,iy-1,iz).u
          +lab(ix+1,iy,iz).v - lab(ix-1,iy,iz).v);
      const Real D13 = .5*factor*(lab(ix,iy,iz+1).u - lab(ix,iy,iz-1).u
          +lab(ix+1,iy,iz).w - lab(ix-1,iy,iz).w);
      const Real D23 = .5*factor*(lab(ix,iy+1,iz).w - lab(ix,iy-1,iz).w
          +lab(ix,iy,iz+1).v - lab(ix,iy,iz-1).v);

      const Real SijTerm =  ( D11*D11 + D22*D22 + D33*D33 + 2*(D12*D12 + D13*D13 + D23*D23) ); // need to multiply this by 2*nu*h^3

      const Real dPdx = factor*(lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p);
      const Real dPdy = factor*(lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p);
      const Real dPdz = factor*(lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p);

      const Real pressTerm = -(dPdx*uGrid + dPdy*vGrid + dPdz*wGrid); // Need to multiply this by h^3

      const Real lapU_1 = invHsqr * (lab(ix+1,iy,iz).u - 3*2*lab(ix,iy,iz).u + lab(ix-1,iy,iz).u + lab(ix,iy+1,iz).u + lab(ix,iy-1,iz).u + lab(ix,iy,iz+1).u + lab(ix,iy,iz-1).u);
      const Real lapU_2 = invHsqr * (lab(ix+1,iy,iz).v - 3*2*lab(ix,iy,iz).v + lab(ix-1,iy,iz).v + lab(ix,iy+1,iz).v + lab(ix,iy-1,iz).v + lab(ix,iy,iz+1).v + lab(ix,iy,iz-1).v);
      const Real lapU_3 = invHsqr * (lab(ix+1,iy,iz).w - 3*2*lab(ix,iy,iz).w + lab(ix-1,iy,iz).w + lab(ix,iy+1,iz).w + lab(ix,iy-1,iz).w + lab(ix,iy,iz+1).w + lab(ix,iy,iz-1).w);
      // Store nu * u \cdot (nabla^2 u)
      // WARNING: multiply here by 0.5 coz later will multiply by dissipFactor=2*nu*hCube
      const Real laplacianTerm= 0.5 * (uGrid*lapU_1 + vGrid*lapU_2 + wGrid*lapU_3); // need to multiply this by 2*nu*h^3

      o(ix,iy,iz).tmpU = (1.0-chi)*(laplacianTerm + SijTerm); // wake
      o(ix,iy,iz).tmpV = (1.0-chi)*(pressTerm); // wake
    }
  }
};


template <typename Lab>
class CoordinatorComputeDissipation : public GenericCoordinator
{
protected:
  const Real nu;
  const int* const step;
  const double* const time;
  Real oldKE=0.0;
public:

  CoordinatorComputeDissipation(FluidGridMPI * _grid, const Real _NU, const int* const _step, const double* const _t) :
    GenericCoordinator(_grid), nu(_NU), step(_step), time(_t)  { }

  void operator()(const double dt)
  {
    check("dissipation - start");
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).tmpU = 0;
        b(ix,iy,iz).tmpV = 0;
        b(ix,iy,iz).tmpW = 0; //zero fields, going to contain Udef
      }
    }

    const Real h = grid->getH();
    const Real ext[3] = {
        h*grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        h*grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        h*grid->getBlocksPerDimension(2)*FluidBlock::sizeZ
    };

    {
      OperatorDissipation kernel(dt,ext,nu);
      compute(kernel);
    }


    double viscous=0.0, press=0.0;
    #pragma omp parallel for schedule(static) reduction (+ : viscous, press )
    for(size_t i=0; i<vInfo.size(); i++)
    {
      Real localViscous = 0.0;
      Real localPress = 0.0;

      BlockInfo info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        localViscous += b(ix,iy,iz).tmpU;
        localPress += b(ix,iy,iz).tmpV;
      }
      viscous += localViscous;
      press += localPress;
    }

    const Real gridH = vInfo[0].h_gridpoint;
    const Real hCube = pow(gridH,3);
    const Real dissipFactor = 2.0*hCube*nu;

    double localSum[2] = {viscous, press};
    double globalSum[2] = {0.0, 0.0};
    MPI_Allreduce(localSum, globalSum, 2, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

    globalSum[0] *= dissipFactor;
    globalSum[1] *= hCube;

    int rank;
    MPI_Comm_rank(grid->getCartComm(),&rank);
    if(rank==0) {

    std::stringstream &fileDissip = logger.get_stream("wakeDissipation.dat");
    if(*step==0)
      fileDissip<<"step_id time viscousTerm pressureTerm total"<<std::endl;

    fileDissip<< *step <<" "<< *time <<" "<< globalSum[0] <<" "<<
	    globalSum[1] <<" "<< globalSum[0]+globalSum[1] <<std::endl;
    }

    check("dissipation - end");
  }

  std::string getName()
  {
    return "Dissipation";
  }
};
#endif
