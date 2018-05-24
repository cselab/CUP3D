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

class NacaMidlineData : public FishMidlineData
{
 protected:

  inline Real _naca_width(const double s, const Real L, const double t_ratio=0.12)
  {
    if(s<0 or s>L) return 0;
    const Real a = 0.2969;
    const Real b =-0.1260;
    const Real c =-0.3516;
    const Real d = 0.2843;
    const Real e =-0.1015;
    //const Real t = 0.12*L;
    const Real t = t_ratio*L;
    const Real p = s/L;

    /*
    if(s>0.99*L){ // Go linear, otherwise trailing edge is not closed - NACA analytical's fault
      const Real temp = 0.99;
      const Real y1 = 5*t* (a*sqrt(temp) +b*temp +c*temp*temp +d*temp*temp*temp + e*temp*temp*temp*temp);
      const Real dydx = (0-y1)/(L-0.99*L);
      return y1 + dydx * (s - 0.99*L);
    }else{ // NACA analytical
      return 5*t* (a*sqrt(p) +b*p +c*p*p +d*p*p*p + e*p*p*p*p);
    }
    */

    return 5*t* (a*std::sqrt(p) +b*p +c*p*p +d*p*p*p + e*p*p*p*p);
  }

 public:
  NacaMidlineData(const int Nm, const double length, const double dx_ext,
    double zExtent, double HoverL=1) : FishMidlineData(Nm,length,1,0,dx_ext)
  {
    #if defined(BC_PERIODICZ)
      // large enough in z-dir such that we for sure fill entire domain
      for(int i=0;i<Nm;++i) height[i] = zExtent;
    #else
      for(int i=0;i<Nm;++i) height[i] = length*HoverL/2;
    #endif
    for(int i=0;i<Nm;++i) width[i]  = _naca_width(rS[i], length);

    computeMidline(0);

    #if 1
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank!=0) return;
      FILE * f = fopen("fish_profile","w");
      for (int i=0; i<Nm; ++i) printf("%g %g %g %g %g\n",
      rX[i],rY[i],rS[i],width[i],height[i]);
      fflush(f); fclose(f);
      printf("Dumped midline\n");
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

  if(!rank)
    printf("Naca: pos=%3.3f, Apitch=%3.3f, Fpitch=%3.3f,Ppitch=%3.3f, "
    "Mpitch=%3.3f, Frow=%3.3f, Arow=%3.3f\n", position[0], Apitch, Fpitch,
    Ppitch, Mpitch, Fheave, Aheave);

  const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
  const double target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
  const double dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
  const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1.)) + 1;
  if (obstacleID) {
    printf("IF3D_NacaOperator WORKS ONLY FOR SINGLE OBSTACLE!\n");
    MPI_Abort(grid->getCartComm(),0);
  }

  myFish = new NacaMidlineData(Nm, length, dx_extension, ext_Z);
}

void IF3D_NacaOperator::update(const int stepID, const double t, const double dt, const Real* Uinf)
{
    _2Dangle  =  Mpitch +          Apitch * std::cos(2*M_PI*(Fpitch*t+Ppitch));
    quaternion[0] = std::cos(0.5*_2Dangle);
    quaternion[1] = 0;
    quaternion[2] = 0;
    quaternion[3] = std::sin(0.5*_2Dangle);

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
    if(bFixFrameOfRef[2]) position[1] += dt*transVel[1];
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

void IF3D_NacaOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{
  //if (!bCreated)
  IF3D_FishOperator::create(step_id, time, dt, Uinf);
}

void IF3D_NacaOperator::finalize(const int step_id,const double time, const double dt, const Real *Uinf)
{
  //if (!bCreated) {
  //  bCreated = true;
    IF3D_FishOperator::finalize(step_id,time, dt, Uinf);
  //} else characteristic_function();
}

void IF3D_NacaOperator::writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock)
{
  #pragma omp parallel
  {
    PutNacaOnBlocks putfish(myFish, position, quaternion);

    #pragma omp for schedule(dynamic)
    for(int i=0; i<vInfo.size(); i++) {
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
}
