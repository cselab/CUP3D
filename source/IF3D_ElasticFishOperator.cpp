//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_ElasticFishOperator.h"
#include "IF3D_FishLibrary.h"
#include "GenericOperator.h"

#include <HDF5Dumper_MPI.h>

#include <cmath>

class ElasticMidlineData : public FishMidlineData
{
 Real * const xOld;
 Real * const yOld;
 Real * const vxOld;
 Real * const vyOld;
 Real * const theta;
 Real * const cosTheta;
 Real * const thetaOld;
 Real * const cosThetaOld;
 public:
  ElasticMidlineData(const double L, const double _h, double zExtent, double t_ratio, double HoverL=1) :
  FishMidlineData(L, 1, 0, _h), xOld(_alloc(Nm)), yOld(_alloc(Nm)), vxOld(_alloc(Nm)),vyOld(_alloc(Nm)), theta(_alloc(Nm)),cosTheta(_alloc(Nm)), thetaOld(_alloc(Nm)),cosThetaOld(_alloc(Nm))
  {
    #if defined(BC_PERIODICZ)
      // large enough in z-dir such that we for sure fill entire domain
      for(int i=0;i<Nm;++i) height[i] = zExtent;
    #else
      for(int i=0;i<Nm;++i) height[i] = length*HoverL/2;
    #endif
    MidlineShapes::naca_width(t_ratio, length, rS, width, Nm);

    computeMidline(0);

    #if 1
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank!=0) return;
      FILE * f = fopen("fish_profile","w");
      for (int i=0; i<Nm; ++i) fprintf(f, "%g %g %g %g %g\n",
      rX[i],rY[i],rS[i],width[i],height[i]);
      fflush(f); fclose(f);
    #endif
  }

  void computeMidline(const double time) override
  {
    rX[0] = rY[0] = vX[0] = vY[0] = 0;
    for(int i=1; i<Nm; ++i) {
      rY[i] = vX[i] = vY[i] = 0;
      rX[i] = rX[i-1] + std::fabs(rS[i]-rS[i-1]);
    }
    _computeMidlineNormals();
/*
      const std::array<double ,6> curvature_points = {
          0, .15*length, .4*length, .65*length, .9*length, length
      };
      const std::array<double ,6> curvature_values = {
          0.82014/length, 1.46515/length, 2.57136/length,
          3.75425/length, 5.09147/length, 5.70449/length
      };
      curvScheduler.transition(time,0,1,curvature_values,curvature_values);
      // query the schedulers for current values
      curvScheduler.gimmeValues(time, curvature_points, Nm, rS, rC, vC);
      // construct the curvature
      for(int i=0; i<Nm; i++) {
        const double darg = 2.*M_PI;
        const double arg  = 2.*M_PI*(time -rS[i]/length) + M_PI*phaseShift;
        rK[i] =   rC[i]*std::sin(arg);
        vK[i] =   vC[i]*std::sin(arg) + rC[i]*std::cos(arg)*darg;
      }

      // solve frenet to compute midline parameters
      IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
      */
  }
};

IF3D_ElasticFishOperator::IF3D_ElasticFishOperator(FluidGridMPI*g, ArgumentParser&p,
  const Real*const u) : IF3D_FishOperator(g, p, u)
{
  const double thickness = p("-thickness").asDouble(0.12); // (NON DIMENSIONAL)
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  myFish = new ElasticMidlineData(length, vInfo[0].h_gridpoint, ext_Z, thickness);
}

void IF3D_ElasticFishOperator::writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock)
{
  #pragma omp parallel
  {
    PutNacaOnBlocks putfish(myFish, position, quaternion);

    #pragma omp for schedule(dynamic)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      const auto pos = segmentsPerBlock.find(info.blockID);
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;

      if(pos != segmentsPerBlock.end()) {
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
          b(ix,iy,iz).tmpU = 0.; //this will be accessed with plus equal

        assert(obstacleBlocks.find(info.blockID) != obstacleBlocks.end());
        ObstacleBlock*const block = obstacleBlocks.find(info.blockID)->second;
        putfish(info, b, block, pos->second);
      }
    }
  }

  #if 0
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic)
    for (int i = 0; i < (int)vInfo.size(); ++i) {
      BlockInfo info = vInfo[i];
      const auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        //b(ix,iy,iz).chi = pos->second->chi[iz][iy][ix];//b(ix,iy,iz).tmpU;
        b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
        b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
        b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
      }
    }
  }
  DumpHDF5_MPI<FluidGridMPI,StreamerVelocityVector>(*grid, 0, 0, "SFD", "./");
  abort();
  #endif
}

void IF3D_ElasticFishOperator::computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump)
{
  IF3D_ObstacleOperator::computeForces(stepID, time, dt, Uinf, NU, bDump);
  // This obstacle requires forces and torques on the midline. War plan:
  // 0) Fetch
  const int Nm = myFish->Nm;
  Real * const fX = myFish->forceX;
  Real * const fY = myFish->forceY;
  Real * const tZ = myFish->torque;
  const Real*const pX = myFish->rX;
  const Real*const pY = myFish->rY;
  // 1) Reset
  std::fill(fX, fX+Nm, 0.0);
  std::fill(fY, fY+Nm, 0.0);
  std::fill(tZ, tZ+Nm, 0.0);
  // 2) Sum surface forces to the closest midline point using section marker
  for (const auto& pos : obstacleBlocks) {
    const ObstacleBlock* const o = pos.second;
    for(int i=0; i<o->nPoints; i++) {
      const int ss = o->ss[i];
      assert(ss>=0 && ss<Nm);
      fX[ss] += o->fX[i];
      fY[ss] += o->fY[i];
      tZ[ss] += (pX[ss]-o->pX[i])*o->fY[i] + (pY[ss]-o->pY[i])*o->fX[i];
    }
  }
  // 3) all reduce across ranks
  MPI_Allreduce(MPI_IN_PLACE, fX, Nm, MPI_DOUBLE, MPI_SUM, grid->getCartComm());
  MPI_Allreduce(MPI_IN_PLACE, fY, Nm, MPI_DOUBLE, MPI_SUM, grid->getCartComm());
  MPI_Allreduce(MPI_IN_PLACE, tZ, Nm, MPI_DOUBLE, MPI_SUM, grid->getCartComm());
}
