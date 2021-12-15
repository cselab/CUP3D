//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Naca.h"
#include "FishLibrary.h"
#include "FishShapes.h"

#include <Cubism/ArgumentParser.h>
#include <Cubism/HDF5Dumper_MPI.h>

#include <cmath>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class NacaMidlineData : public FishMidlineData
{
 Real * const rK;
 Real * const vK;
 Real * const rC;
 Real * const vC;
 public:
  NacaMidlineData(const double L, const double _h, double zExtent,
    double t_ratio, double HoverL=1) : FishMidlineData(L, 1, 0, _h),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm))
  {
    for(int i=0;i<Nm;++i) height[i] = length*HoverL/2;
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
  rZ[i] = 0;
  vZ[i] = 0;
      }
      #pragma omp parallel for schedule(static)
      for(int i=0; i<Nm-1; i++) {
        const double ds = rS[i+1]-rS[i];
        const double tX = rX[i+1]-rX[i];
        const double tY = rY[i+1]-rY[i];
        const double tVX = vX[i+1]-vX[i];
        const double tVY = vY[i+1]-vY[i];
        norX[i] = -tY/ds;
        norY[i] =  tX/ds;
        norZ[i] =  0.0;
        vNorX[i] = -tVY/ds;
        vNorY[i] =  tVX/ds;
        vNorZ[i] = 0.0;
        binX[i] =  0.0;
        binY[i] =  0.0;
        binZ[i] =  1.0;
        vBinX[i] = 0.0;
        vBinY[i] = 0.0;
        vBinZ[i] = 0.0;
      }
      norX[Nm-1] = norX[Nm-2];
      norY[Nm-1] = norY[Nm-2];
      norZ[Nm-1] = norZ[Nm-2];
      vNorX[Nm-1] = vNorX[Nm-2];
      vNorY[Nm-1] = vNorY[Nm-2];
      vNorZ[Nm-1] = vNorZ[Nm-2];
      binX[Nm-1] = binX[Nm-2];
      binY[Nm-1] = binY[Nm-2];
      binZ[Nm-1] = binZ[Nm-2];
      vBinX[Nm-1] = vBinX[Nm-2];
      vBinY[Nm-1] = vBinY[Nm-2];
      vBinZ[Nm-1] = vBinZ[Nm-2];
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

Naca::Naca(SimulationData&s, ArgumentParser&p) : Fish(s, p)
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

  if(!sim.rank)
    printf("Naca: pos=%3.3f, Apitch=%3.3f, Fpitch=%3.3f,Ppitch=%3.3f, "
    "Mpitch=%3.3f, Frow=%3.3f, Arow=%3.3f\n", position[0], Apitch, Fpitch,
    Ppitch, Mpitch, Fheave, Aheave);
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  myFish = new NacaMidlineData(length, sim.hmin, sim.extent[2], thickness);
}

void Naca::update()
{
  _2Dangle =  Mpitch + Apitch * std::cos( 2*M_PI * (Fpitch*sim.time + Ppitch) );
  quaternion[0] = std::cos(0.5*_2Dangle);
  quaternion[1] = 0;
  quaternion[2] = 0;
  quaternion[3] = std::sin(0.5*_2Dangle);

  absPos[0] += sim.dt * transVel[0];
  absPos[1] += sim.dt * transVel[1];
  absPos[2] += sim.dt * transVel[2];

  position[0] += sim.dt * ( transVel[0] + sim.uinf[0] );
  // if user wants to keep airfoil in the mid plane then we just integrate
  // relative velocity (should be 0), otherwise we know that y velocity
  // is sinusoidal, therefore we can just use analytical form
  if(bFixFrameOfRef[1])
    position[1] += sim.dt * ( transVel[1] + sim.uinf[1] );
  else
    position[1] = sim.extent[1]/2 + Aheave * std::cos(2*M_PI*Fheave*sim.time);
  position[2] += sim.dt * ( transVel[2] + sim.uinf[2] );

  _writeComputedVelToFile();
}

void Naca::computeVelocities()
{
  Obstacle::computeVelocities();

  // x velocity can be either fixed from the start (then we follow the obst op
  // pattern) or self propelled, here we do not touch it.
  const Real argv = 2*M_PI * Fheave * sim.time;
  transVel[1] = -2*M_PI * Fheave * Aheave * std::sin(argv);
  transVel[2] = 0;
  angVel[0] = 0;
  angVel[1] = 0;
  const Real arga = 2*M_PI * ( Fpitch * sim.time + Ppitch );
  angVel[2] = -2*M_PI * Fpitch * Apitch * std::sin(arga);
}

using intersect_t = std::vector<std::vector<VolumeSegment_OBB*>>;
void Naca::writeSDFOnBlocks(std::vector<VolumeSegment_OBB> & vSegments)
{
  #pragma omp parallel
  {
    PutNacaOnBlocks putfish(myFish, position, quaternion);
    #pragma omp for
    for (size_t j=0 ; j < MyBlockIDs.size(); j++)
    {
      std::vector<VolumeSegment_OBB*> S;
      for (size_t k = 0 ; k < MySegments[j].size() ; k++)
      {
        VolumeSegment_OBB*const ptr  = & vSegments[MySegments[j][k]];
        S.push_back(ptr);
      }
      if(S.size() > 0)
      {
        ObstacleBlock*const block = obstacleBlocks[MyBlockIDs[j].blockID];
        putfish(MyBlockIDs[j].h,
                MyBlockIDs[j].origin_x,
                MyBlockIDs[j].origin_y,
                MyBlockIDs[j].origin_z, block, S);
      }
    }
  }
}

CubismUP_3D_NAMESPACE_END
