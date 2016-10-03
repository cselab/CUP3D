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

        void _prepareRotation2D(const Real angle)
        {
        	Rmatrix2D[0][0] = Rmatrix2D[1][1] = std::cos(angle);
        	Rmatrix2D[0][1] = -std::sin(angle);
        	Rmatrix2D[1][0] = -Rmatrix2D[0][1];
        }

        inline void _subtractAngularVelocity(const Real angvel, const Real x, const Real y, Real & vx, Real & vy) const
        {
        	vx += angvel*y;
        	vy -= angvel*x;
        }

        void _computeWidthsHeights()
        {
            for(int i=0;i<Nm;++i) {
                width[i]  = Fish::width(rS[i],length);
                height[i] = Fish::height(rS[i],length);
            }
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

        void _computeMidlineNormals()
        {
#pragma omp parallel for
            for(int i=0; i<Nm-1; i++) {
                const Real ds = rS[i+1]-rS[i];
                const Real tX = rX[i+1]-rX[i];
                const Real tY = rY[i+1]-rY[i];
                const Real tVX = vX[i+1]-vX[i];
                const Real tVY = vY[i+1]-vY[i];
                norX[i] = -tY/ds;
                norY[i] =  tX/ds;
                vNorX[i] = -tVY/ds;
                vNorY[i] =  tVX/ds;
            }
            norX[Nm-1] = norX[Nm-2];
            norY[Nm-1] = norY[Nm-2];
            vNorX[Nm-1] = vNorX[Nm-2];
            vNorY[Nm-1] = vNorY[Nm-2];
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

public:
        FishMidlineData(const int Nm, const Real len, const Real Tp, const Real phase, const Real dx_ext):
        	Nm(Nm),length(len),Tperiod(Tp),phaseShift(phase),rS(_alloc(Nm)),rX(_alloc(Nm)),rY(_alloc(Nm)),
			vX(_alloc(Nm)),vY(_alloc(Nm)),norX(_alloc(Nm)),norY(_alloc(Nm)),vNorX(_alloc(Nm)),vNorY(_alloc(Nm)),
			width(_alloc(Nm)),height(_alloc(Nm)),iFishStart(4*NPPEXT),iFishEnd(Nm-1-4*NPPEXT)
        {
        	// extension_info contains number of extension points and extension dx
            const int Nextension = 4*NPPEXT; // up to 3dx on each side (to get proper interpolation up to 2dx)
        	const int Next = Nextension; // number of points per extension
        	const int Nint = Nm -2*Next; // number of interior points

        	// extension head
        	for(int i=0;i<Next;++i)
        		rS[i] = 0.0 - (Next- i) * dx_ext;
        	// interior points
        	for(int i=0;i<Nint;++i)
        		rS[i+Next] = length * 0.5 * (1.0 - std::cos(i * M_PI/((Real)Nint-1))); // cosine: more points near head and tail
        	// rS[i] = i*length/((Real)Nint-1); // linear: equally distributed points
        	// extension tail
        	for(int i=0;i<Next;++i)
        		rS[i+Nint+Next] = length + (i + 1)*dx_ext;
        	_computeWidthsHeights();
        }

        ~FishMidlineData()
        {
            _dealloc(rS);
            _dealloc(rX);
            _dealloc(rY);
            _dealloc(vX);
            _dealloc(vY);
            _dealloc(norX);
            _dealloc(norY);
            _dealloc(vNorX);
            _dealloc(vNorY);
            _dealloc(height);
            _dealloc(width);
        }

        Real integrateLinearMomentum(Real CoM[2], Real vCoM[2])
        {   // already worked out the integrals for r, theta on paper
            // remaining integral done with composite trapezoidal rule
            // minimize rhs evaluations --> do first and last point separately
            Real _vol(0), _cmx(0), _cmy(0), _lmx(0), _lmy(0);
#pragma omp parallel for reduction(+:_vol,_cmx,_cmy,_lmx,_lmy)
            for(int i=0;i<Nm;++i) {
                const Real ds = (i==0) ? rS[1]-rS[0] :
                		((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
                const Real fac1 = _integrationFac1(i);
                const Real fac2 = _integrationFac2(i);
                _vol += 0.5*fac1*ds;
                _cmx += 0.5*(rX[i]*fac1 + norX[i]*fac2)*ds;
                _cmy += 0.5*(rY[i]*fac1 + norY[i]*fac2)*ds;
                _lmx += 0.5*(vX[i]*fac1 + vNorX[i]*fac2)*ds;
                _lmy += 0.5*(vY[i]*fac1 + vNorY[i]*fac2)*ds;
            }

            vol=_vol*M_PI;
            CoM[0]=_cmx*M_PI;
            CoM[1]=_cmy*M_PI;
            linMom[0]=_lmx*M_PI;
            linMom[1]=_lmy*M_PI;

            assert(vol> std::numeric_limits<Real>::epsilon());
            const Real ivol = 1.0/vol;

            CoM[0]*=ivol;
            CoM[1]*=ivol;
            vCoM[0]=linMom[0]*ivol;
            vCoM[1]=linMom[1]*ivol;
            //printf("%f %f %f %f %f\n",CoM[0],CoM[1],vCoM[0],vCoM[1], vol);
            return vol;
        }

        void integrateAngularMomentum(Real & angVel)
        {
            // assume we have already translated CoM and vCoM to nullify linear momentum

            // already worked out the integrals for r, theta on paper
            // remaining integral done with composite trapezoidal rule
            // minimize rhs evaluations --> do first and last point separately
            Real _J(0), _am(0);

            #pragma omp parallel for reduction(+:_J,_am)
            for(int i=0;i<Nm;++i) {
            	const Real ds = (i==0) ? rS[1]-rS[0] :
            			((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
                const Real fac1 = _integrationFac1(i);
                const Real fac2 = _integrationFac2(i);
                const Real fac3 = _integrationFac3(i);
                double tmp_J, tmp_M;
                tmp_M  = (rX[i]*vY[i] - rY[i]*vX[i])*fac1;
                tmp_M += (rX[i]*vNorY[i] - rY[i]*vNorX[i] + vY[i]*norX[i] - vX[i]*norY[i])*fac2;
                tmp_M += (norX[i]*vNorY[i] - norY[i]*vNorX[i])*fac3;
                _am += 0.5*tmp_M*ds;
                tmp_J  = (rX[i]*rX[i] + rY[i]*rY[i])*fac1;
                tmp_J += 2.0*(rX[i]*norX[i] + rY[i]*norY[i])*fac2;
                //tmpSum += (norX[idx]*norX[idx]+norY[idx]*norY[idx])*fac3;
                tmp_J += fac3;
                _J += 0.5*tmp_J*ds;
            }

            J=_J*M_PI;
            angMom=_am*M_PI;
            assert(J>std::numeric_limits<Real>::epsilon());
            angVel = angMom/J;
        }

        void changeToCoMFrameLinear(const Real CoM_internal[2], const Real vCoM_internal[2])
        {
            for(int i=0;i<Nm;++i) {
                rX[i]-=CoM_internal[0];
                rY[i]-=CoM_internal[1];
                vX[i]-=vCoM_internal[0];
                vY[i]-=vCoM_internal[1];
            }
        }

        void changeToCoMFrameAngular(const Real theta_internal, const Real angvel_internal)
        {
            _prepareRotation2D(theta_internal);
#pragma omp parallel for
            for(int i=0;i<Nm;++i) {
                _rotate2D(rX[i],rY[i]);
                _rotate2D(vX[i],vY[i]);
                _subtractAngularVelocity(angvel_internal, rX[i], rY[i], vX[i], vY[i]);
            }
            _computeMidlineNormals();
        }

        virtual void computeMidline(const Real time) = 0;
        virtual void _correctTrajectory(const Real dtheta, const Real time, const Real dt) {}
        virtual void execute(const Real time, const Real l_tnext, const vector<Real>& input) {}
    };

    struct CurvatureDefinedFishData : FishMidlineData
    {
    protected:
    	Schedulers::ParameterSchedulerVector<6> curvScheduler;
        Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
    	Schedulers::ParameterSchedulerVector<6> adjustScheduler;
    	Real * const rK;
    	Real * const vK;
    	Real * const rC;
    	Real * const vC;
    	Real * const rB;
    	Real * const vB;
    	Real * const rA;
    	Real * const vA;
    	Real l_Tp, time0, timeshift;

    public:

    	CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
    	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), l_Tp(Tperiod), timeshift(0), time0(0),
		rK(_alloc(Nm)),vK(_alloc(Nm)),rC(_alloc(Nm)),vC(_alloc(Nm)),rA(_alloc(Nm)),vA(_alloc(Nm)),rB(_alloc(Nm)),vB(_alloc(Nm))
    	{    	}

    	void _correctTrajectory(const Real dtheta, const Real time, const Real dt) override
    	{
    		std::array<Real,6> tmp_curv = std::array<Real,6>();
    		for (int i=0; i<tmp_curv.size(); ++i) {tmp_curv[i] = dtheta/M_PI;}
    		adjustScheduler.transition(time,time,time+2*dt,tmp_curv);
    	}

        void execute(const Real time, const Real l_tnext, const vector<Real>& input) override
        {
            if (input.size()>1) {
                baseScheduler.Turn(input[0], l_tnext);
                //first, shift time to  previous turn node
                timeshift += (l_tnext-time0)/l_Tp;
                time0 = l_tnext;
                l_Tp = Tperiod*(1.+input[1]);
            } else if (input.size()>0) {
                baseScheduler.Turn(input[0], l_tnext);
            }
        }

    	~CurvatureDefinedFishData()
    	{
    		_dealloc(rK);
    		_dealloc(vK);
    		_dealloc(rC);
    		_dealloc(vC);
    		_dealloc(rB);
    		_dealloc(vB);
    		_dealloc(rA);
    		_dealloc(vA);
    	}

    	void computeMidline(const Real time)
    	{
    		const Real _1oL = 1./length;
    		const std::array<Real ,6> curvature_values = {
    				0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
					3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
    		};
    		const std::array<Real ,6> curvature_points = {
    			0., .15*length, .4*length, .65*length, .9*length, length
    		};
    		const std::array<Real ,7> baseline_points = {
    				1.00, 0.75, 0.50, 0.25, 0.00, -0.25, -0.50
    		};
    		const std::array<Real, 6> curvature_zeros = std::array<Real, 6>();
    		curvScheduler.transition(time,0.0,Tperiod,curvature_zeros,curvature_values);

    		// query the schedulers for current values
    		curvScheduler.gimmeValues(time, curvature_points, Nm, rS, rC, vC);
    		baseScheduler.gimmeValues(time, l_Tp, length, baseline_points, Nm, rS, rB, vB);
    		adjustScheduler.gimmeValues(time, curvature_points, Nm, rS, rA, vA);

    		// construct the curvature
    		const Real _1oT = 1./l_Tp;
			for(unsigned int i=0; i<Nm; i++) {
				const Real darg = 2.*M_PI* _1oT;
				const Real arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL) + M_PI*phaseShift;
				rK[i] = rC[i]*(std::sin(arg) + rB[i] + rA[i]);
				vK[i] = vC[i]*(std::sin(arg) + rB[i] + rA[i]) + rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]);
			}

#if 1==0
    		{ // we dump the profile points
    			FILE * f = fopen("stefan.dat","a");
    			std::array<Real, 6> curv,base;
    			curvScheduler.ParameterScheduler<6>::gimmeValues(time, curv);
    			baseScheduler.ParameterScheduler<6>::gimmeValues(time, base);
    			fprintf(f,"%10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",time,curv[0],curv[1],curv[2],curv[3],curv[4],curv[5],base[0],base[1],base[2],base[3],base[4],base[5]);
    			fclose(f);
    		}
    		{ // we dump the profile
    			FILE * f = fopen("stefan_profile","w");
    			for(int i=0;i<Nm;++i)
    				fprintf(f,"%d %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i],vNorX[i],vNorY[i],width[i]);
    			fclose(f);
    		}
#endif

    		// solve frenet to compute midline parameters
    		IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
    	}
    };

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

        void _computeMidlineCoordinates(const Real time)
        {
        	const Real rampFac = rampFactorSine(time, Tperiod);
        	rX[0] = 0.0;
        	rY[0] =   rampFac*midline(rS[0], time, length, Tperiod, phaseShift);

        	for(int i=1;i<Nm;++i) {
        		rY[i]=rampFac*midline(rS[i], time, length, Tperiod, phaseShift);
        		const Real dy = rY[i]-rY[i-1];
        		const Real ds = rS[i] - rS[i-1];
        		const Real dx = std::sqrt(ds*ds-dy*dy);
        		rX[i] = rX[i-1] + dx;
        	}
        }

        void _computeMidlineVelocities(const Real time)
        {
        	const Real rampFac =    rampFactorSine(time, Tperiod);
        	const Real rampFacVel = rampFactorVelSine(time, Tperiod);

        	vX[0] = 0.0; //rX[0] is constant
        	vY[0] = rampFac*midlineVel(rS[0],time,length,Tperiod, phaseShift) +
        			rampFacVel*midline(rS[0],time,length,Tperiod, phaseShift);

        	for(int i=1;i<Nm;++i) {
        		vY[i]=rampFac*midlineVel(rS[i],time,length,Tperiod, phaseShift) +
        			  rampFacVel*midline(rS[i],time,length,Tperiod, phaseShift);
        		const Real dy = rY[i]-rY[i-1];
        		const Real dx = rX[i]-rX[i-1];
        		const Real dVy = vY[i]-vY[i-1];
        		assert(dx>0); // has to be, otherwise y(s) is multiple valued for a given s
        		vX[i] = vX[i-1] - dy/dx * dVy; // use ds^2 = dx^2 + dy^2 --> ddx = -dy/dx*ddy
        	}
        }

    public:
    	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
    	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext)
    	{
#ifdef BBURST
    		ifstream reader("burst_coast_carling_params.txt");
    		if (reader.is_open()) {
    			reader >> t0;
    			reader >> t1;
    			reader >> t2;
    			reader >> t3;
    			reader.close();
    		} else {
    			cout << "Could not open params.txt" << endl;
    			abort();
    		}
#endif
    	}

        void computeMidline(const Real time)
        {
            _computeMidlineCoordinates(time);
            _computeMidlineVelocities(time);
            _computeMidlineNormals();
#ifndef NDEBUG
    		// we dump the profile
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			if (rank!=0) return;
			FILE * f = fopen("fish_profile","w");
			for(int i=0;i<Nm;++i)
				fprintf(f,"%d %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e %10.10e\n",
						i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i],vNorX[i],vNorY[i],width[i],height[i]);
			fclose(f);
			printf("Dumped midline\n");
#endif
        }
    };
}

class IF3D_CarlingFishOperator: public IF3D_ObstacleOperator
{
	Fish::CarlingFishMidlineData * myFish;
    Real Tperiod, phaseShift, phase, sim_time, sim_dt;
    
    Real volume_internal, J_internal, CoM_internal[2], vCoM_internal[2], theta_internal, angvel_internal, angvel_internal_prev, CoM_interpolated[3];
    bool randomStart, randomActions, bSpiral, useLoadedActions, bCorrectTrajectory;
    //mt19937 * gen; //TODO... what to do? shared seed?
    //smarties:
    Real Tstartlearn, GoalDX, new_curv, old_curv, new_Tp, adjTh;
    vector<vector<Real>> loadedActions;
    int  nActions;

    //Real signLastTurn;
    //bool bObstacleBlocksFilled;
    //int nTurnActions, nPauseActions;
    //const bool bFixToPlanar; // restrict translation to xy plane, and rotation to z-axis

public:
	
    IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
    : IF3D_ObstacleOperator(grid, parser), theta_internal(0.0), angvel_internal(0.0), sim_time(0.0), sim_dt(0.0), adjTh(adjTh), myFish(nullptr)
	{
        volume=0;
        for(int i=0;i<3;i++) transVel[i]=0;
        for(int i=0;i<3;i++) angVel[i]=0;
        for(int i=0;i<6;i++) J[i]=0;
        _parseArguments(parser);
        const Real target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
        const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
        const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG)+1;
        printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
        fflush(0);
        // multiple of NPPSEG: TODO why?
        myFish = new Fish::CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension);
    }
    
    ~IF3D_CarlingFishOperator()
    {
    	if(myFish not_eq nullptr) delete myFish;
    }

	void save(const int step_id, const Real t, std::string filename = std::string()) override;
	void restart(const Real t, std::string filename = std::string()) override;
    void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    void getCenterOfMass(Real CM[3]) const override;
    void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
