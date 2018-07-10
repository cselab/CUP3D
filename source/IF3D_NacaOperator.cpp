//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_NacaOperator.h"
#include "IF3D_FishLibrary.h"
#include "GenericOperator.h"

#include <HDF5Dumper_MPI.h>

#include <cmath>

class NacaMidlineData : public FishMidlineData
{
 Real * const rK;
 Real * const vK;
 Real * const rC;
 Real * const vC;
 public:
  NacaMidlineData(const double L, const double _h, double zExtent, double t_ratio, double HoverL=1) :
  FishMidlineData(L, 1, 0, _h),rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm))
  {
    #if defined(BC_PERIODICZ)
      // large enough in z-dir such that we for sure fill entire domain
      for(int i=0;i<Nm;++i) height[i] = zExtent;
    #else
      for(int i=0;i<Nm;++i) height[i] = length*HoverL/2;
    #endif
    MidlineShapes::naca_width(t_ratio, length, rS, width, Nm);

    computeMidline(0.0, 0.0);

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

  void computeMidline(const double time, const double dt) override
  {
    #if 1
      rX[0] = rY[0] = vX[0] = vY[0] = 0;
      for(int i=1; i<Nm; ++i) {
        rY[i] = vX[i] = vY[i] = 0;
        rX[i] = rX[i-1] + std::fabs(rS[i]-rS[i-1]);
      }
      _computeMidlineNormals();
    #else // 2d stefan swimmer
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
    #endif
  }
};

IF3D_NacaOperator::IF3D_NacaOperator(FluidGridMPI*g, ArgumentParser&p,
  const Real*const u) : IF3D_FishOperator(g, p, u), bCreated(false)
{
  absPos[0] = 0;
  #if 1
      Apitch = p("-Apitch").asDouble(0.0); //aplitude of sinusoidal pitch angle
      Fpitch = p("-Fpitch").asDouble(0.0); //frequency
      Ppitch = p("-Ppitch").asDouble(0.0); //phase wrt to rowing motion
      Mpitch = p("-Mpitch").asDouble(0.0); //mean angle
      Fheave = p("-Fheave").asDouble(0.0); //frequency of rowing motion
      Aheave = p("-Aheave").asDouble(0.0); //amplitude (NON DIMENSIONAL)
  #else
      ifstream reader("params.txt");
      if (reader.is_open()) {
        Apitch=0.0; Fpitch=0.0; Ppitch=0.0; Mpitch=0.0; Fheave=0.0; Aheave=0.0;
        reader >> Apitch; //aplitude of sinusoidal pitch angle
        reader >> Fpitch; //frequency
        reader >> Ppitch; //phase wrt to rowing motion
        reader >> Mpitch; //mean angle
        reader >> Fheave; //frequency of rowing motion
        reader >> Aheave; //amplitude of rowing motion
        reader.close();
      } else {
        cout << "Could not open params.txt" << endl;
        abort();
      }
  #endif
  Aheave *= length;
  isSelfPropelled = Aheave > 0;

  const double thickness = p("-thickness").asDouble(0.12); // (NON DIMENSIONAL)

  if(!rank)
    printf("Naca: pos=%3.3f, Apitch=%3.3f, Fpitch=%3.3f,Ppitch=%3.3f, "
    "Mpitch=%3.3f, Frow=%3.3f, Arow=%3.3f\n", position[0], Apitch, Fpitch,
    Ppitch, Mpitch, Fheave, Aheave);
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  myFish = new NacaMidlineData(length, vInfo[0].h_gridpoint, ext_Z, thickness);
}

void IF3D_NacaOperator::update(const int stepID, const double t, const double dt, const Real* Uinf)
{
    _2Dangle  =  Mpitch +          Apitch * std::cos(2*M_PI*(Fpitch*t+Ppitch));
    quaternion[0] = std::cos(0.5*_2Dangle);
    quaternion[1] = 0;
    quaternion[2] = 0;
    quaternion[3] = std::sin(0.5*_2Dangle);

    #if 0
      _2Dangle = 0.6;
      quaternion[0] = std::cos(0.5*_2Dangle);
      quaternion[1] = 0;
      quaternion[2] = 0;
      quaternion[3] = std::sin(0.5*_2Dangle);
      const double one = std::sqrt( quaternion[0]*quaternion[0]
                                   +quaternion[1]*quaternion[1]
                                   +quaternion[2]*quaternion[2]
                                   +quaternion[3]*quaternion[3]);
      quaternion[0] /= one+2e-16;
      quaternion[1] /= one+2e-16;
      quaternion[2] /= one+2e-16;
      quaternion[3] /= one+2e-16;
    #endif

    Real velx_tot = transVel[0]-Uinf[0];
    Real vely_tot = transVel[1]-Uinf[1];
    Real velz_tot = transVel[2]-Uinf[2];
    absPos[0] += dt*velx_tot;
    absPos[1] += dt*vely_tot;
    absPos[2] += dt*velz_tot;

    position[0] += dt*transVel[0];
    // if user wants to keep airfoil in the mid plane then we just integrate
    // relative velocity (should be 0), otherwise we know that y velocity
    // is sinusoidal, therefore we can just use analytical form
    if(bFixFrameOfRef[1]) position[1] += dt*transVel[1];
    else position[1] = ext_Y/2 + Aheave * std::cos(2*M_PI*Fheave*t);
    position[2] += dt*transVel[2];

    sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                      _2Dangle, velx_tot, vely_tot, angVel[2]);
    tOld = t;
    _writeComputedVelToFile(stepID, t, Uinf);
}

void IF3D_NacaOperator::computeVelocities(const Real* Uinf)
{
  IF3D_ObstacleOperator::computeVelocities(Uinf);

  // x velocity can be either fixed from the start (then we follow the obst op
  // pattern) or self propelled, here we do not touch it.
  transVel[1] = -2*M_PI*Fheave*Aheave *std::sin(2*M_PI*Fheave*tOld);
  transVel[2] = 0;
  angVel[0] = 0;
  angVel[1] = 0;
  angVel[2] = -2*M_PI*Fpitch*Apitch *std::sin(2*M_PI*(Fpitch*tOld+Ppitch));
}

void IF3D_NacaOperator::writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock)
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
