//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "../SimulationData.h"
#include "../operators/Operator.h"
#include "PoissonSolverBase.h"
#include "PoissonSolverAMRKernels.h"
#include <Cubism/BlockInfo.h>
#include <vector>
#include <cassert>
#include <cstring>

namespace cubismup3d {

class ComputeLHS : public Operator
{
  struct KernelLHSPoisson
  {
    const SimulationData & sim;
    KernelLHSPoisson(const SimulationData&s) : sim(s) {}
    const std::vector<BlockInfo> & lhsInfo = sim.lhsInfo();
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;
    const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};

    void operator()(const ScalarLab & lab, const BlockInfo& info) const
    {
      ScalarBlock & __restrict__ o = (*sim.lhs)(info.blockID);
      const Real h = info.h; 
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        #ifdef PRESERVE_SYMMETRY
        o(x,y,z).s = h * ( ConsistentSum(lab(x-1,y,z).s + lab(x+1,y,z).s, lab(x,y-1,z).s + lab(x,y+1,z).s , lab(x,y,z-1).s + lab(x,y,z+1).s) - 6.0*lab(x,y,z).s );
        #else
        o(x,y,z) = h*( lab(x-1,y,z) + lab(x+1,y,z) + lab(x,y-1,z) + lab(x,y+1,z) +lab(x,y,z-1) + lab(x,y,z+1) - 6.0*lab(x,y,z));
        #endif
      }

      BlockCase<ScalarBlock> * tempCase = (BlockCase<ScalarBlock> *)(lhsInfo[info.blockID].auxiliary);

      if (tempCase == nullptr) return; //no flux corrections needed for this block

      ScalarElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
      ScalarElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
      ScalarElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
      ScalarElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
      ScalarElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
      ScalarElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;

      if (faceXm != nullptr)
      {
        const int x = 0;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXm[y + Ny*z] = h*(lab(x,y,z) - lab(x-1,y,z));
      }
      if (faceXp != nullptr)
      {
        const int x = Nx-1;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXp[y + Ny*z] = h*(lab(x,y,z) - lab(x+1,y,z));
      }
      if (faceYm != nullptr)
      {
        const int y = 0;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYm[x + Nx*z] = h*(lab(x,y,z) - lab(x,y-1,z));
      }
      if (faceYp != nullptr)
      {
        const int y = Ny-1;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYp[x + Nx*z] = h*(lab(x,y,z) - lab(x,y+1,z));
      }
      if (faceZm != nullptr)
      {
        const int z = 0;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZm[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z-1));
      }
      if (faceZp != nullptr)
      {
        const int z = Nz-1;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZp[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z+1));
      }

    }
  };

  public:
  ComputeLHS(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    Real avgP = 0;
    int index = -1;
    MPI_Request request;

    const std::vector<BlockInfo> & vInfo_lhs = sim.lhsInfo();
    const std::vector<BlockInfo> & vInfo_z   = sim.presInfo();
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;

    if (sim.bMeanConstraint <= 2 && sim.bMeanConstraint > 0)
    {
      #pragma omp parallel for reduction(+ : avgP)
      for(size_t i=0; i<vInfo_z.size(); ++i)
      {
        const ScalarBlock & __restrict__ Z  = (*sim.pres)(i);
        const Real h3 = vInfo_z[i].h*vInfo_z[i].h*vInfo_z[i].h;
        if (vInfo_z[i].index[0] == 0 && vInfo_z[i].index[1] == 0 &&  vInfo_z[i].index[2] == 0) index = i;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          avgP += Z(x,y,z).s*h3;
      }
      MPI_Iallreduce(MPI_IN_PLACE, &avgP, 1, MPI_Real, MPI_SUM, sim.comm, &request);
    }

    compute<ScalarLab>(KernelLHSPoisson(sim),sim.pres,sim.lhs);

    if (sim.bMeanConstraint == 0) return;

    if (sim.bMeanConstraint <= 2 && sim.bMeanConstraint > 0)
    {
      MPI_Waitall(1,&request,MPI_STATUSES_IGNORE);
      if (sim.bMeanConstraint == 1 && index != -1)
      {
        ScalarBlock & __restrict__ LHS  = (*sim.lhs)(index);
        LHS(0,0,0).s = avgP;
      }
      else if (sim.bMeanConstraint == 2)
      {
        #pragma omp parallel for
        for(size_t i=0; i<vInfo_lhs.size(); ++i)
	{
          ScalarBlock & __restrict__ LHS = (*sim.lhs)(i);
          const Real h3 = vInfo_lhs[i].h*vInfo_lhs[i].h*vInfo_lhs[i].h;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
            LHS(x,y,z).s += avgP*h3;
	}
      }
    }
    else // > 2
    {
      #pragma omp parallel for
      for(size_t i=0; i<vInfo_lhs.size(); ++i)
      {
        ScalarBlock & __restrict__ LHS = (*sim.lhs)(i);
        const ScalarBlock & __restrict__ Z = (*sim.pres)(i);
        if (vInfo_lhs[i].index[0] == 0 && vInfo_lhs[i].index[1] == 0 && vInfo_lhs[i].index[2] == 0) LHS(0,0,0).s = Z(0,0,0).s;
      }
    }
  }
  std::string getName() { return "ComputeLHS"; }
};

class PoissonSolverAMR : public PoissonSolverBase
{
 protected:
  SimulationData& sim;
  ComputeLHS findLHS;

#if 1
  void _preconditioner(const std::vector<Real> & input, std::vector<Real> & output)
  {
    auto &  zInfo         = sim.pres->getBlocksInfo(); //used for preconditioning
    const size_t Nblocks  = zInfo.size();
    const int BSX         = VectorBlock::sizeX;
    const int BSY         = VectorBlock::sizeY;
    const int BSZ         = VectorBlock::sizeZ;

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ bb = (*sim.pres)(i);
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        bb(ix,iy,iz).s = input[j];
      }
    }

    #pragma omp parallel
    {
      cubismup3d::poisson_kernels::getZImplParallel(sim.presInfo());
    }

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      const ScalarBlock & __restrict__ bb = (*sim.pres)(i);
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        output[j] = bb(ix,iy,iz).s;;
      }
    }

  }

  void _lhs(std::vector<Real> & input, std::vector<Real> & output)
  {
    auto &  zInfo         = sim.pres->getBlocksInfo(); //used for preconditioning
    auto & AxInfo         = sim.lhs ->getBlocksInfo(); //will store the LHS result
    const size_t Nblocks  = zInfo.size();
    const int BSX         = VectorBlock::sizeX;
    const int BSY         = VectorBlock::sizeY;
    const int BSZ         = VectorBlock::sizeZ;

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ zz = *(ScalarBlock*) zInfo[i].ptrBlock;
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        zz(ix,iy,iz).s  = input[j];
      }
    }

    findLHS(0);

    #pragma omp parallel for
    for (size_t i=0; i < Nblocks; i++)
    {
      ScalarBlock & __restrict__ Ax = *(ScalarBlock*) AxInfo[i].ptrBlock;
      for(int iz=0; iz<BSZ; iz++)
      for(int iy=0; iy<BSY; iy++)
      for(int ix=0; ix<BSX; ix++)
      {
        const int j = i*BSX*BSY*BSZ+iz*BSX*BSY+iy*BSX+ix;
        output[j]   = Ax(ix,iy,iz).s;
      }
    }
  }

  std::vector<Real> b;
  std::vector<Real> phat;
  std::vector<Real> rhat;
  std::vector<Real> shat;
  std::vector<Real> what;
  std::vector<Real> zhat;
  std::vector<Real> qhat;
  std::vector<Real> s;
  std::vector<Real> w;
  std::vector<Real> z;
  std::vector<Real> t;
  std::vector<Real> v;
  std::vector<Real> q;
  std::vector<Real> r;
  std::vector<Real> y;
  std::vector<Real> x;
  std::vector<Real> r0;
  std::vector<Real> x_opt;

#else //non-pipelined version

  void getZ()
  {
    #pragma omp parallel
    {
      cubismup3d::poisson_kernels::getZImplParallel(sim.presInfo());
    }
  }
  std::vector<Real> x   ;
  std::vector<Real> r   ;
  std::vector<Real> p   ;
  std::vector<Real> v   ;
  std::vector<Real> s   ;
  std::vector<Real> rhat;

  size_t _dest(const BlockInfo &info , const int iz, const int iy, const int ix) const
  {
    return BlockType::sizeX * ( BlockType::sizeY * (info.blockID * BlockType::sizeZ  + iz) + iy) + ix;
  }

#endif

 public:
  PoissonSolverAMR(SimulationData& ss): sim(ss),findLHS(ss){}
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 
  void solve();
};

}//namespace cubismup3d






