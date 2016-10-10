//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_CarlingFishOperator__
#define __IncompressibleFluids3D__IF3D_CarlingFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

namespace Fish
{
    struct CarlingFishMidlineData : FishMidlineData
    {
    protected:
        //burst-coast:
        Real t0, t1, t2, t3;

    	inline Real rampFactorSine(const Real t, const Real T) const
    	{
    		return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
    	}

    	inline Real rampFactorVelSine(const Real t, const Real T) const
    	{
    		return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
    	}

    	inline Real midline(const Real s, const Real t, const Real L, const Real T, const Real phaseShift) const
    	{
    		const Real arg = 2.0*M_PI*(s/L - t/T + phaseShift);
        	const Real fac = 0.1212121212121212;
        	const Real inv = 0.03125;

#ifdef BBURST
    		Real f;
    		Real tcoast = TSTART;
    		if (t>=TSTART) {
    			const Real bct = t0 + t1 + t2 + t3;
    			const Real phase = std::floor((t-TSTART)/bct);
    			tcoast = TSTART + phase*bct;
    		}
    		Real tfreeze = tcoast + t0;
    		Real tburst = tfreeze + t1;
    		Real tswim = tburst + t2;

    		if (t<tcoast) {
    			f = 1.0;
    		} else if (t<tfreeze) {
    			const Real d = (t-tcoast)/(tfreeze-tcoast);
    			f = 1 - 3*d*d + 2*d*d*d;
    		} else if (t<tburst) {
    			f = 0.0;
    		} else if (t<tswim) {
    			const Real d = (t-tburst)/(tswim-tburst);
    			f = 3*d*d - 2*d*d*d;
    		} else {
    			f = 1.0;
    		}
    		return f * fac * (s + inv*L)*std::sin(arg);
#else
    		return fac * (s + inv*L)*std::sin(arg);
#endif
    	}

    	inline Real midlineVel(const Real s, const Real t, const Real L, const Real T, const Real phaseShift) const
    	{
        	const Real arg = 2.0*M_PI*(s/L - t/T + phaseShift);
        	const Real fac = 0.1212121212121212;
        	const Real inv = 0.03125;

#ifdef BBURST
        	Real f,df;
        	Real tcoast = TSTART;
        	if (t>=TSTART) {
        		const Real bct = t0 + t1 + t2 + t3;
        		const Real phase = std::floor((t-TSTART)/bct);
        		tcoast = TSTART + phase*bct;
        	}
        	Real tfreeze = tcoast + t0;
        	Real tburst = tfreeze + t1;
        	Real tswim = tburst + t2;

        	if (t<tcoast) {
        		f = 1.0;
        		df = 0.0;
        	} else if (t<tfreeze) {
        		const Real d = (t-tcoast)/(tfreeze-tcoast);
        		f = 1 - 3*d*d + 2*d*d*d;
        		df = 6*(d*d - d)/(tfreeze-tcoast);
        	} else if (t<tburst) {
        		f = 0.0;
        		df = 0.0;
        	} else if (t<tswim) {
        		const Real d = (t-tburst)/(tswim-tburst);
        		f = 3*d*d - 2*d*d*d;
        		df = 6*(d - d*d)/(tswim-tburst);
        	} else {
        		f = 1.0;
        		df = 0.0;
        	}

        	return fac*(s + inv*L)*(df*std::sin(arg) - f*(2.0*M_PI/T)*std::cos(arg));
#else
    		return - fac*(s + inv*L)*(2.0*M_PI/T)*std::cos(arg);
#endif //Burst-coast
    	}

        void _computeMidlineCoordinates(const Real time);
        void _computeMidlineVelocities(const Real time);

    public:
    	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext);
        void computeMidline(const Real time);
    };
}

class IF3D_CarlingFishOperator: public IF3D_FishOperator
{
public:
	
    IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser);
    ~IF3D_CarlingFishOperator();
	void save(const int step_id, const Real t, std::string filename = std::string()) override;
	void restart(const Real t, std::string filename = std::string()) override;
    void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    void getCenterOfMass(Real CM[3]) const override;
    void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
