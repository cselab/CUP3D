//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_VortexOperator__
#define __IncompressibleFluids3D__IF3D_VortexOperator__

#include <cmath>
#include <array>
#include "IF3D_ObstacleOperator.h"

class IF3D_VortexOperator: public IF3D_ObstacleOperator
{
protected:
    double v_max;
    bool created;

public:
    IF3D_VortexOperator(FluidGridMPI* g, ArgumentParser& p, const Real*const u);
    virtual void update(const int step_id, const double t, const double dt, const Real *Uinf) override;
    virtual void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
    void computeVelocities(const Real* Uinf) override;
    void computeForces(const int stepID, const double time, const double dt,
                       const Real* Uinf, const double NU, const bool bDump) override;
  void Accept(ObstacleVisitor * visitor) override;
  void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void execute(const int iAgent, const double time, const vector<double>a) override;
  void interpolateOnSkin(const double time, const int stepID, bool dumpWake=false) override;
  void getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts) override;
  void computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_Vortex__) */
