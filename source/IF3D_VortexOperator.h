//
//  IF3D_VortexOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_VortexOperator__
#define __IncompressibleFluids3D__IF3D_VortexOperator__

#include <cmath>
#include <array>
#include "IF3D_ObstacleOperator.h"

class IF3D_VortexOperator: public IF3D_ObstacleOperator
{
protected:
    Real v_max;
		bool created;

public:
    IF3D_VortexOperator(FluidGridMPI * grid, ArgumentParser & parser);
    virtual void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    virtual void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
		void computeVelocities(const Real* Uinf) override;
		void computeForces(const int stepID, const Real time, const Real dt,
                       const Real* Uinf, const Real NU, const bool bDump) override;
	void Accept(ObstacleVisitor * visitor) override;
	void finalize(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
	void execute(Communicator * comm, const int iAgent, const Real time, const int iLabel) override;
	void interpolateOnSkin(const Real time, const int stepID, bool dumpWake=false) override;
	void getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts) override;
	void computeDiagnostics(const int stepID, const Real time, const Real* Uinf, const Real lambda) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_Vortex__) */
