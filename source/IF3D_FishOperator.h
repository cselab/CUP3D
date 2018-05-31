//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_FishOperator__
#define __IncompressibleFluids3D__IF3D_FishOperator__

#include "IF3D_ObstacleOperator.h"

class FishMidlineData;
struct VolumeSegment_OBB;
typedef std::map<int, std::vector<VolumeSegment_OBB>> mapBlock2Segs;
typedef std::vector<VolumeSegment_OBB> aryVolSeg;

class IF3D_FishOperator: public IF3D_ObstacleOperator
{
protected:
  FishMidlineData * myFish = nullptr;
  //phaseShift=0, phase=0,
  double Tperiod=0, sim_time=0, sim_dt=0;
  double volume_internal=0, J_internal=0, CoM_internal[2]={0,0}, vCoM_internal[2]={0,0};
  double theta_internal=0, angvel_internal=0, angvel_internal_prev=0;
  double CoM_interpolated[3] = {0,0,0}, angvel_integral[3] = {0,0,0};
  double adjTh=0, adjDy=0, followX=0, followY=0;
  bool bCorrectTrajectory=false;
  //const Real* ptrUinf_copy = nullptr;

  void integrateMidline();
  virtual void writeSDFOnBlocks(const mapBlock2Segs& segmentsPerBlock);
  void apply_pid_corrections(const double time, const double dt, const Real *Uinf);
  aryVolSeg prepare_vSegments();
  //override to create special obstacle blocks for local force balances:
  virtual mapBlock2Segs prepare_segPerBlock(const aryVolSeg&vSegments);

public:
  IF3D_FishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  ~IF3D_FishOperator();
  void save(const int step_id, const double t, std::string filename = std::string()) override;
  void restart(const double t, std::string filename = std::string()) override;

  virtual void update(const int step_id, const double t, const double dt, const Real *Uinf) override;

  void getCenterOfMass(double CM[3]) const override;
  void interpolateOnSkin(const double time, const int stepID, bool dumpWake=false) override;

  virtual void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
  virtual void computeChi(const int step_id, const double time, const double dt, const Real *Uinf, int& mpi_status) override;
  virtual void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;

  void getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL, Real*& pYL,
    Real*& pXU, Real*& pYU, int& Npts) override;

  //  void computeVelocities(const Real Uinf[3]) override
  //  {
  //    computeVelocities_forced(Uinf);
  //  }
  // void setTranslationVelocity(double UT[3]) override  { }
};


#endif /* defined(__IncompressibleFluids3D__IF3D_Fish__) */
