//
//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"
#ifdef _ACCFFT_
#include "PoissonSolverScalarFFTW_ACC.h"
typedef PoissonSolverScalarFFTW_ACC<FluidGridMPI, StreamerDiv> PressureSolver;
#else
#ifdef _UNBOUNDED_FFT_
#include "PoissonSolverScalarFFTW_cyclicConvolution.h"
#else
#include "PoissonSolverScalarFFTW.h"
#endif /* _UNBOUNDED_FFT_ */
//#ifdef FFTW_FFT
typedef PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> PressureSolver;
//#else
//typedef PoissonSolverScalarFFTW_DCT_MPI<FluidGridMPI,StreamerDiv>PressureSolver;
//#endif
#endif


struct PressureObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * grid;
  vector<BlockInfo> vInfo;

  PressureObstacleVisitor(FluidGridMPI * _grid) : grid(_grid)
  {
    vInfo = grid->getBlocksInfo();
  }

   void visit(IF3D_ObstacleOperator* const obstacle)
   {
     #pragma omp parallel
     {
       const std::map<int,ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
       #pragma omp for schedule(dynamic)
       for (int i = 0; i < (int)vInfo.size(); ++i) {
         BlockInfo info = vInfo[i];
         const auto pos = obstblocks.find(info.blockID);
         if(pos == obstblocks.end()) continue;

         FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

         for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
         for(int iy=0; iy<FluidBlock::sizeY; ++iy)
         for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
           // what if multiple obstacles share a block??
           // let's plus equal and wake up during the night to stress about it
           b(ix,iy,iz).tmpU += pos->second->udef[iz][iy][ix][0];
           b(ix,iy,iz).tmpV += pos->second->udef[iz][iy][ix][1];
           b(ix,iy,iz).tmpW += pos->second->udef[iz][iy][ix][2];
         }
       }
     }
   }
};

class OperatorDivergenceMinusDivTmpU : public GenericLabOperator
{
 private:
  const double dt;
  PressureSolver*const solver;

 public:
  OperatorDivergenceMinusDivTmpU(double _dt, PressureSolver*const s)
  : dt(_dt), solver(s)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 6, 1,2,3,5,6,7);
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
  }
  ~OperatorDivergenceMinusDivTmpU() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const size_t offset = solver->_offset_ext(info);
    const Real fac = 0.5/info.h_gridpoint/dt;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      // Poisson solver reads field p for the rhs
      const Real dU    = lab(ix+1,iy  ,iz  ).u - lab(ix-1,iy  ,iz  ).u;
      const Real dV    = lab(ix  ,iy+1,iz  ).v - lab(ix  ,iy-1,iz  ).v;
      const Real dW    = lab(ix  ,iy  ,iz+1).w - lab(ix  ,iy  ,iz-1).w;
      const Real dUdef = lab(ix+1,iy  ,iz  ).tmpU - lab(ix-1,iy  ,iz  ).tmpU;
      const Real dVdef = lab(ix  ,iy+1,iz  ).tmpV - lab(ix  ,iy-1,iz  ).tmpV;
      const Real dWdef = lab(ix  ,iy  ,iz+1).tmpW - lab(ix  ,iy  ,iz-1).tmpW;
      const Real ret = fac*(dU+dV+dW -lab(ix,iy,iz).chi*(dUdef+dVdef+dWdef));
      // const Real ret = fac*(dU+dV+dW);  // Use this to kill internal divergence.
      //o(ix,iy,iz).p = ret;
      solver->_cub2fftw(offset, iz, iy, ix, ret);
    }
  }
};

class OperatorDivergenceMinusDivTmpU2ndOrder : public GenericLabOperator
{
 private:
  double dt;
  PressureSolver*const solver;

 public:
  OperatorDivergenceMinusDivTmpU2ndOrder(double _dt, PressureSolver*const s) : dt(_dt), solver(s)
  {
    stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 6, 1,2,3,5,6,7);
    stencil_start[0] = -2; stencil_start[1] = -2;  stencil_start[2] = -2;
    stencil_end[0] = 3;  stencil_end[1] = 3;  stencil_end[2] = 3;
  }
  ~OperatorDivergenceMinusDivTmpU2ndOrder() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const size_t offset = solver->_offset_ext(info);
    const Real factor = 1./(12*info.h_gridpoint*dt);

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      // Poisson solver reads field p for the rhs
      const Real dU     = lab(ix+1,iy  ,iz  ).u    - lab(ix-1,iy  ,iz  ).u;
      const Real dV     = lab(ix  ,iy+1,iz  ).v    - lab(ix  ,iy-1,iz  ).v;
      const Real dW     = lab(ix  ,iy  ,iz+1).w    - lab(ix  ,iy  ,iz-1).w;
      const Real dU2    = lab(ix+2,iy  ,iz  ).u    - lab(ix-2,iy  ,iz  ).u;
      const Real dV2    = lab(ix  ,iy+2,iz  ).v    - lab(ix  ,iy-2,iz  ).v;
      const Real dW2    = lab(ix  ,iy  ,iz+2).w    - lab(ix  ,iy  ,iz-2).w;
      const Real dUdef  = lab(ix+1,iy  ,iz  ).tmpU - lab(ix-1,iy  ,iz  ).tmpU;
      const Real dVdef  = lab(ix  ,iy+1,iz  ).tmpV - lab(ix  ,iy-1,iz  ).tmpV;
      const Real dWdef  = lab(ix  ,iy  ,iz+1).tmpW - lab(ix  ,iy  ,iz-1).tmpW;
      const Real dUdef2 = lab(ix+2,iy  ,iz  ).tmpU - lab(ix-2,iy  ,iz  ).tmpU;
      const Real dVdef2 = lab(ix  ,iy+2,iz  ).tmpV - lab(ix  ,iy-2,iz  ).tmpV;
      const Real dWdef2 = lab(ix  ,iy  ,iz+2).tmpW - lab(ix  ,iy  ,iz-2).tmpW;

      const Real tmp1 = 8*(dU + dV + dW) - (dU2 + dV2 + dW2);
      const Real tmp2 = 8*(dUdef + dVdef + dWdef) - (dUdef2 + dVdef2 + dWdef2);
      const Real ret  = factor*(tmp1 -o(ix,iy,iz).chi*tmp2);
      solver->_cub2fftw(offset, iz, iy, ix, ret);
      //o(ix,iy,iz).p = ret;
    }
  }
};

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
    const Real fac = - 0.5 * dt/ info.h_gridpoint;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      // p contains the pressure correction after the Poisson solver
      o(ix,iy,iz).u += fac * (lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p);
      o(ix,iy,iz).v += fac * (lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p);
      o(ix,iy,iz).w += fac * (lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p);
    }
  }
};

class OperatorGradP2ndOrder : public GenericLabOperator
{
 private:
  const double dt;
  const Real extent[3];

 public:
  OperatorGradP2ndOrder(double _dt,const Real ext[3]) :
  dt(_dt),extent{ext[0],ext[1],ext[2]}
  {
    stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 1, 4);
    stencil_start[0] = -2; stencil_start[1] = -2; stencil_start[2] = -2;
    stencil_end[0] = 3; stencil_end[1] = 3; stencil_end[2] = 3;
  }

  ~OperatorGradP2ndOrder() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = -dt / info.h_gridpoint / 12;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const Real dPx = lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p;
      const Real dPy = lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p;
      const Real dPz = lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p;
      const Real dPx2= lab(ix+2,iy,iz).p - lab(ix-2,iy,iz).p;
      const Real dPy2= lab(ix,iy+2,iz).p - lab(ix,iy-2,iz).p;
      const Real dPz2= lab(ix,iy,iz+2).p - lab(ix,iy,iz-2).p;
      // p contains the pressure correction after the Poisson solver
      o(ix,iy,iz).u += fac * (8*dPx - dPx2);
      o(ix,iy,iz).v += fac * (8*dPy - dPy2);
      o(ix,iy,iz).w += fac * (8*dPz - dPz2);
    }
  }
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
 protected:
  PressureSolver pressureSolver;
  IF3D_ObstacleVector** const obstacleVector;

 public:
  CoordinatorPressure(FluidGridMPI * _grid, IF3D_ObstacleVector** const myobstacles) :
    GenericCoordinator(_grid),
    pressureSolver(*_grid),
    obstacleVector(myobstacles)
  { }

  void operator()(const double dt)
  {
    check("pressure - start");
    const int nthreads = omp_get_max_threads();

    #pragma omp parallel
    {
      const int N = vInfo.size();
      #pragma omp for schedule(static)
      for(int i=0; i<N; i++) {
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
    }

    const Real h = grid->getH();
    const Real ext[3] = {
        h*grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        h*grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        h*grid->getBlocksPerDimension(2)*FluidBlock::sizeZ
    };
     //store deformation velocities onto tmp fields
    ObstacleVisitor * pressureVisitor = new PressureObstacleVisitor(grid);
    (*obstacleVector)->Accept(pressureVisitor); // accept you son of a french cow
    delete pressureVisitor;

    {   //place onto p: ( div u^(t+1) - div u^* ) / dt
      //where i want div u^(t+1) to be equal to div udef
      vector<OperatorDivergenceMinusDivTmpU*> diff(nthreads, nullptr);
      for(int i=0;i<nthreads;++i)
        diff[i] = new OperatorDivergenceMinusDivTmpU(dt, &pressureSolver);

      compute<OperatorDivergenceMinusDivTmpU>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    pressureSolver.solve();

    { //pressure correction dudt* = - grad P / rho
      vector<OperatorGradP*> diff(nthreads, nullptr);
      for(int i=0;i<nthreads;++i) diff[i] = new OperatorGradP(dt, ext);

      compute<OperatorGradP>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    check("pressure - end");
  }

  string getName()
  {
    return "Pressure";
  }
};
#endif
