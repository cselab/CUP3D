//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_FishOperator__
#define __IncompressibleFluids3D__IF3D_FishOperator__

#include <cmath>
#include <array>
#include "IF2D_Frenet.h"
#include "IF3D_ObstacleOperator.h"

const int NPPSEG = 50.; //was 100
const int NPPEXT = 3; //was 3
const int TGTPPB = 4.; //was 2 i think
const int TSTART = 2.;

namespace Fish
{
	inline Real width(const Real s, const Real L)
	{
		if(s<0 or s>L) return 0;
		const Real sb = .04*L;
		const Real st = .95*L;
		const Real wt = .01*L;
		const Real wh = .04*L;

		return (s<sb ? std::sqrt(2.0*wh*s-s*s) :
				(s<st ? wh-(wh-wt)*std::pow((s-sb)/(st-sb),2) :
						(wt * (L-s)/(L-st))));
	}

	inline Real height(const Real s, const Real L)
	{
		if(s<0 or s>L) return 0;
		const Real a=0.51*L;
		const Real b=0.08*L;
		return b*std::sqrt(1 - std::pow((s-a)/a,2));
	}

    struct FishMidlineData
    {
    public:
        const int Nm;
        Real * const rS; // arclength discretization points
        Real * const rX; // coordinates of midline discretization points
        Real * const rY;
        Real * const vX; // midline discretization velocities
        Real * const vY;
        Real * const norX; // normal vector to the midline discretization points
        Real * const norY;
        Real * const vNorX;
        Real * const vNorY;
        Real * const width;
        Real * const height;
        Real linMom[2], vol, J, angMom; // for diagnostics
        // start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
        const int iFishStart, iFishEnd;

    protected:
        const Real length;
        const Real Tperiod;
        const Real phaseShift;
        Real Rmatrix2D[2][2];
        Real Rmatrix3D[3][3];

        inline void _rotate2D(Real &x, Real &y) const
        {
        	const Real p[2] = {x,y};
        	x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
        	y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
        }

        inline void _translateAndRotate2D(const Real pos[2], Real &x, Real &y) const
        {
        	const Real p[2] = {
        			x-pos[0],
					y-pos[1]
        	};
        	// rotate
        	x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
        	y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
        }

        inline Real _d_ds(const int idx, const Real * const vals, const int maxidx) const
        {
            if(idx==0)
                return (vals[idx+1]-vals[idx])/(rS[idx+1]-rS[idx]);
            else if(idx==maxidx-1)
                return (vals[idx]-vals[idx-1])/(rS[idx]-rS[idx-1]);
            else
                return 0.5*( (vals[idx+1]-vals[idx])/(rS[idx+1]-rS[idx]) + (vals[idx]-vals[idx-1])/(rS[idx] - rS[idx-1]) );
        }

        Real * _alloc(const int N)
        {
            return new Real[N];
        }

        void _dealloc(Real * ptr)
        {
            if(ptr not_eq nullptr) {
                delete [] ptr;
                ptr=nullptr;
            }
        }

        inline Real _integrationFac1(const int idx) const
        {
            return width[idx]*height[idx];
        }

        inline Real _integrationFac2(const int idx) const
        {
            const Real dnorXi = _d_ds(idx, norX, Nm);
            const Real dnorYi = _d_ds(idx, norY, Nm);
            return 0.25*std::pow(width[idx],3)*height[idx]*(dnorXi*norY[idx] - dnorYi*norX[idx]);
        }

        inline Real _integrationFac3(const int idx) const
        {
            const Real drXi = _d_ds(idx, rX, Nm);
            const Real drYi = _d_ds(idx, rY, Nm);
            // return 0.25*std::pow(width[idx],3)*height[idx]*(drXi*norY[idx] - drYi*norX[idx]);
            return 0.25*std::pow(width[idx],3)*height[idx];
        }

        void _prepareRotation2D(const Real angle);
        void _computeWidthsHeights();
        void _computeMidlineNormals();

public:
        FishMidlineData(const int Nm, const Real len, const Real Tp, const Real phase, const Real dx_ext);
        ~FishMidlineData();
        Real integrateLinearMomentum(Real CoM[2], Real vCoM[2]);
        void integrateAngularMomentum(Real & angVel);
        void changeToCoMFrameLinear(const Real CoM_internal[2], const Real vCoM_internal[2]);
        void changeToCoMFrameAngular(const Real theta_internal, const Real angvel_internal);
        virtual void computeMidline(const Real time) = 0;
        virtual void _correctTrajectory(const Real dtheta, const Real time, const Real dt) {}
        virtual void execute(const Real time, const Real l_tnext, const vector<Real>& input) {}
    };
}

class IF3D_FishOperator: public IF3D_ObstacleOperator
{
protected:
	Fish::FishMidlineData * myFish;
    Real Tperiod, phaseShift, phase, sim_time, sim_dt;
    Real volume_internal, J_internal, CoM_internal[2], vCoM_internal[2];
    Real theta_internal, angvel_internal, angvel_internal_prev, CoM_interpolated[3];
    Real Tstartlearn, GoalDX, new_curv, old_curv, new_Tp, adjTh, angvel_integral[3];
    bool randomStart, bCorrectTrajectory;
    int  nActions;

public:
	
    IF3D_FishOperator(FluidGridMPI * grid, ArgumentParser & parser);
    ~IF3D_FishOperator();
	void save(const int step_id, const Real t, std::string filename = std::string()) override;
	void restart(const Real t, std::string filename = std::string()) override;
    void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    void getCenterOfMass(Real CM[3]) const override;
    void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_Fish__) */
