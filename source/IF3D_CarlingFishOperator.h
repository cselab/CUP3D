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
	inline double width(const double s, const Real L)
	{
		if(s<0 or s>L) return 0;
		const double sb = .04*L;
		const double st = .95*L;
		const double wt = .01*L;
		const double wh = .04*L;

		return (s<sb ? std::sqrt(2.0*wh*s-s*s) :
				(s<st ? wh-(wh-wt)*std::pow((s-sb)/(st-sb),2) :
						(wt * (L-s)/(L-st))));
	}

	inline double height(const double s, const Real L)
	{
		if(s<0 or s>L) return 0;
		const double a=0.51*L;
		const double b=0.08*L;
		return b*std::sqrt(1 - std::pow((s-a)/a,2));
	}

    struct FishMidlineData
    {
    public:
        const int Nm;
        double * const rS; // arclength discretization points
        double * const rX; // coordinates of midline discretization points
        double * const rY;
        double * const vX; // midline discretization velocities
        double * const vY;
        double * const norX; // normal vector to the midline discretization points
        double * const norY;
        double * const vNorX;
        double * const vNorY;
        double * const width;
        double * const height;
        double linMom[2], vol, J, angMom; // for diagnostics
        // start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
        const int iFishStart, iFishEnd;

    protected:
        const Real length;
        const Real Tperiod;
        const Real phaseShift;
        double Rmatrix2D[2][2];
        double Rmatrix3D[3][3];

        inline void _rotate2D(double &x, double &y) const
        {
        	const double p[2] = {x,y};
        	x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
        	y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
        }

        inline void _translateAndRotate2D(const double pos[2], double &x, double &y) const
        {
        	const double p[2] = {
        			x-pos[0],
					y-pos[1]
        	};
        	// rotate
        	x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
        	y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
        }

        void _prepareRotation2D(const double angle)
        {
        	Rmatrix2D[0][0] = Rmatrix2D[1][1] = std::cos(angle);
        	Rmatrix2D[0][1] = -std::sin(angle);
        	Rmatrix2D[1][0] = -Rmatrix2D[0][1];
        }

        inline void _subtractAngularVelocity(const double angvel, const double x, const double y, double & vx, double & vy) const
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

        inline double _d_ds(const int idx, const double * const vals, const int maxidx) const
        {
            if(idx==0)
                return (vals[idx+1]-vals[idx])/(rS[idx+1]-rS[idx]);
            else if(idx==maxidx-1)
                return (vals[idx]-vals[idx-1])/(rS[idx]-rS[idx-1]);
            else
                return 0.5*( (vals[idx+1]-vals[idx])/(rS[idx+1]-rS[idx]) + (vals[idx]-vals[idx-1])/(rS[idx] - rS[idx-1]) );
        }

        double * _alloc(const int N)
        {
            return new double[N];
        }

        void _dealloc(double * ptr)
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
                const double ds = rS[i+1]-rS[i];
                const double tX = rX[i+1]-rX[i];
                const double tY = rY[i+1]-rY[i];
                const double tVX = vX[i+1]-vX[i];
                const double tVY = vY[i+1]-vY[i];
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

        inline double _integrationFac1(const int idx) const
        {
            return width[idx]*height[idx];
        }

        inline double _integrationFac2(const int idx) const
        {
            const double dnorXi = _d_ds(idx, norX, Nm);
            const double dnorYi = _d_ds(idx, norY, Nm);
            return 0.25*std::pow(width[idx],3)*height[idx]*(dnorXi*norY[idx] - dnorYi*norX[idx]);
        }

        inline double _integrationFac3(const int idx) const
        {
            const double drXi = _d_ds(idx, rX, Nm);
            const double drYi = _d_ds(idx, rY, Nm);
            // return 0.25*std::pow(width[idx],3)*height[idx]*(drXi*norY[idx] - drYi*norX[idx]);
            return 0.25*std::pow(width[idx],3)*height[idx];
        }

        inline void _updateVolumeIntegration(const double fac1, const double ds, double & vol) const
        {
            vol+=0.5*fac1*ds;
        }

        inline void _updateCoMIntegration(const int idx, const double fac1, const double fac2, const double ds, double CoM[2]) const
        {
            CoM[0] += 0.5*(rX[idx]*fac1 + norX[idx]*fac2)*ds;
            CoM[1] += 0.5*(rY[idx]*fac1 + norY[idx]*fac2)*ds;
        }

        inline void _updateLinMomIntegration(const int idx, const double fac1, const double fac2, const double ds, double linMom[2]) const
        {
            linMom[0] += 0.5*(vX[idx]*fac1 + vNorX[idx]*fac2)*ds;
            linMom[1] += 0.5*(vY[idx]*fac1 + vNorY[idx]*fac2)*ds;
        }

        inline void _updateAngMomIntegration(const int idx, const double fac1, const double fac2, const double fac3, const double ds, double & angMom) const
        {
            double tmpSum = 0.0;
            tmpSum += (rX[idx]*vY[idx] - rY[idx]*vX[idx])*fac1;
            tmpSum += (rX[idx]*vNorY[idx] - rY[idx]*vNorX[idx] + vY[idx]*norX[idx] - vX[idx]*norY[idx])*fac2;
            tmpSum += (norX[idx]*vNorY[idx] - norY[idx]*vNorX[idx])*fac3;
            angMom += 0.5*tmpSum*ds;
        }

        inline void _updateJIntegration(const int idx, const double fac1, const double fac2, const double fac3, const double ds, double & J) const
        {
            double tmpSum = 0.0;
            tmpSum += (rX[idx]*rX[idx] + rY[idx]*rY[idx])*fac1;
            tmpSum += 2.0*(rX[idx]*norX[idx] + rY[idx]*norY[idx])*fac2;
            // tmpSum += (norX[idx]*norX[idx]+norY[idx]*norY[idx])*fac3;
            tmpSum += fac3;
            J += 0.5*tmpSum*ds;
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
        		rS[i+Next] = length * 0.5 * (1.0 - std::cos(i * M_PI/((double)Nint-1))); // cosine: more points near head and tail
        	// rS[i] = i*length/((double)Nint-1); // linear: equally distributed points
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

        double integrateLinearMomentum(double CoM[2], double vCoM[2])
        {
            CoM[0]=0.0;
            CoM[1]=0.0;
            vol=0;
            linMom[0]=linMom[1]=0;

            // already worked out the integrals for r, theta on paper
            // remaining integral done with composite trapezoidal rule
            // minimize rhs evaluations --> do first and last point separately
            {
                const double ds0 = rS[1]-rS[0];
                const double fac1 = _integrationFac1(0);
                const double fac2 = _integrationFac2(0);
                _updateVolumeIntegration(fac1,ds0,vol);
                _updateCoMIntegration(0, fac1,fac2,ds0,CoM);
                _updateLinMomIntegration(0, fac1,fac2,ds0,linMom);
            }

            for(int i=1;i<Nm-1;++i)
            {
                const double ds = rS[i+1]-rS[i-1];
                const double fac1 = _integrationFac1(i);
                const double fac2 = _integrationFac2(i);
                _updateVolumeIntegration(fac1,ds,vol);
                _updateCoMIntegration(i, fac1,fac2,ds,CoM);
                _updateLinMomIntegration(i, fac1,fac2,ds,linMom);
            }

            {
                const double dse = rS[Nm-1]-rS[Nm-2];
                const double fac1 = _integrationFac1(Nm-1);
                const double fac2 = _integrationFac2(Nm-1);
                _updateVolumeIntegration(fac1,dse,vol);
                _updateCoMIntegration(Nm-1, fac1,fac2,dse,CoM);
                _updateLinMomIntegration(Nm-1, fac1,fac2,dse,linMom);
            }

            vol*=M_PI;
            CoM[0]*=M_PI;
            CoM[1]*=M_PI;
            linMom[0]*=M_PI;
            linMom[1]*=M_PI;

            assert(vol> std::numeric_limits<double>::epsilon());
            const double ivol = 1.0/vol;

            CoM[0]*=ivol;
            CoM[1]*=ivol;
            vCoM[0]=linMom[0]*ivol;
            vCoM[1]=linMom[1]*ivol;
            //printf("%f %f %f %f %f\n",CoM[0],CoM[1],vCoM[0],vCoM[1], vol);
            return vol;
        }

        void integrateAngularMomentum(double & angVel)
        {
            // assume we have already translated CoM and vCoM to nullify linear momentum
            J=0;
            angMom=0;

            // already worked out the integrals for r, theta on paper
            // remaining integral done with composite trapezoidal rule
            // minimize rhs evaluations --> do first and last point separately
            {
                const double ds0 = rS[1]-rS[0];
                const double fac1 = _integrationFac1(0);
                const double fac2 = _integrationFac2(0);
                const double fac3 = _integrationFac3(0);
                _updateJIntegration(0, fac1, fac2, fac3, ds0, J);
                _updateAngMomIntegration(0, fac1, fac2, fac3, ds0, angMom);
            }

            for(int i=1;i<Nm-1;++i)
            {
                const double ds = rS[i+1]-rS[i-1];
                const double fac1 = _integrationFac1(i);
                const double fac2 = _integrationFac2(i);
                const double fac3 = _integrationFac3(i);
                _updateJIntegration(i, fac1, fac2, fac3, ds, J);
                _updateAngMomIntegration(i, fac1, fac2, fac3, ds, angMom);
            }

            {
                const double dse = rS[Nm-1]-rS[Nm-2];
                const double fac1 = _integrationFac1(Nm-1);
                const double fac2 = _integrationFac2(Nm-1);
                const double fac3 = _integrationFac3(Nm-1);
                _updateJIntegration(Nm-1, fac1, fac2, fac3, dse, J);
                _updateAngMomIntegration(Nm-1, fac1, fac2, fac3, dse, angMom);
            }

            J*=M_PI;
            angMom*=M_PI;
            assert(J>std::numeric_limits<double>::epsilon());
            angVel = angMom/J;
        }

        void changeToCoMFrameLinear(const double CoM_internal[2], const double vCoM_internal[2])
        {
            for(int i=0;i<Nm;++i) {
                rX[i]-=CoM_internal[0];
                rY[i]-=CoM_internal[1];
                vX[i]-=vCoM_internal[0];
                vY[i]-=vCoM_internal[1];
            }
        }

        void changeToCoMFrameAngular(const double theta_internal, const double angvel_internal)
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

        virtual void computeMidline(const double time) = 0;
        virtual void _correctTrajectory(const double dtheta, const double time, const double dt) {}
        virtual void execute(const double time, const double l_tnext, const vector<Real>& input) {}
    };

    struct CurvatureDefinedFishData : FishMidlineData
    {
    protected:
    	Schedulers::ParameterSchedulerVector<6> curvScheduler;
        Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
    	Schedulers::ParameterSchedulerVector<6> adjustScheduler;
    	double * const rK;
    	double * const vK;
    	double * const rC;
    	double * const vC;
    	double * const rB;
    	double * const vB;
    	double * const rA;
    	double * const vA;
    	double l_Tp, time0, timeshift;

    public:

    	CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
    	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), l_Tp(Tperiod), timeshift(0), time0(0),
		rK(_alloc(Nm)),vK(_alloc(Nm)),rC(_alloc(Nm)),vC(_alloc(Nm)),rA(_alloc(Nm)),vA(_alloc(Nm)),rB(_alloc(Nm)),vB(_alloc(Nm))
    	{    	}

    	void _correctTrajectory(const double dtheta, const double time, const double dt) override
    	{
    		std::array<double,6> tmp_curv = std::array<double,6>();
    		for (int i=0; i<tmp_curv.size(); ++i) {tmp_curv[i] = dtheta/M_PI;}
    		adjustScheduler.transition(time,time,time+2*dt,tmp_curv);
    	}

        void execute(const double time, const double l_tnext, const vector<Real>& input) override
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

    	void computeMidline(const double time)
    	{
    		const double _1oL = 1./length;
    		const std::array<double ,6> curvature_values = {
    				0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
					3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
    		};
    		const std::array<double ,6> curvature_points = {
    			0., .15*length, .4*length, .65*length, .9*length, length
    		};
    		const std::array<double ,7> baseline_points = {
    				1.00, 0.75, 0.50, 0.25, 0.00, -0.25, -0.50
    		};
    		const std::array<double, 6> curvature_zeros = std::array<double, 6>();
    		curvScheduler.transition(time,0.0,Tperiod,curvature_zeros,curvature_values);

    		// query the schedulers for current values
    		curvScheduler.gimmeValues(time, curvature_points, Nm, rS, rC, vC);
    		baseScheduler.gimmeValues(time, l_Tp, length, baseline_points, Nm, rS, rB, vB);
    		adjustScheduler.gimmeValues(time, curvature_points, Nm, rS, rA, vA);

    		// construct the curvature
    		const double _1oT = 1./l_Tp;
			for(unsigned int i=0; i<Nm; i++) {
				const double darg = 2.*M_PI* _1oT;
				const double arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL) + M_PI*phaseShift;
				rK[i] = rC[i]*(std::sin(arg) + rB[i] + rA[i]);
				vK[i] = vC[i]*(std::sin(arg) + rB[i] + rA[i]) + rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]);
			}

#if 1==0
    		{ // we dump the profile points
    			FILE * f = fopen("stefan.dat","a");
    			std::array<double, 6> curv,base;
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
        double t0, t1, t2, t3;

    	inline double rampFactorSine(const Real t, const Real T) const
    	{
    		return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
    	}

    	inline double rampFactorVelSine(const Real t, const Real T) const
    	{
    		return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
    	}

    	inline double midline(const double s, const Real t, const Real L, const Real T, const Real phaseShift) const
    	{
    		const double arg = 2.0*M_PI*(s/L - t/T + phaseShift);
        	const double fac = 0.1212121212121212;
        	const double inv = 0.03125;

#ifdef BBURST
    		double f;
    		double tcoast = TSTART;
    		if (t>=TSTART) {
    			const double bct = t0 + t1 + t2 + t3;
    			const double phase = std::floor((t-TSTART)/bct);
    			tcoast = TSTART + phase*bct;
    		}
    		double tfreeze = tcoast + t0;
    		double tburst = tfreeze + t1;
    		double tswim = tburst + t2;

    		if (t<tcoast) {
    			f = 1.0;
    		} else if (t<tfreeze) {
    			const double d = (t-tcoast)/(tfreeze-tcoast);
    			f = 1 - 3*d*d + 2*d*d*d;
    		} else if (t<tburst) {
    			f = 0.0;
    		} else if (t<tswim) {
    			const double d = (t-tburst)/(tswim-tburst);
    			f = 3*d*d - 2*d*d*d;
    		} else {
    			f = 1.0;
    		}
    		return f * fac * (s + inv*L)*std::sin(arg);
#else
    		return fac * (s + inv*L)*std::sin(arg);
#endif
    	}

    	inline double midlineVel(const double s, const Real t, const Real L, const Real T, const Real phaseShift) const
    	{
        	const double arg = 2.0*M_PI*(s/L - t/T + phaseShift);
        	const double fac = 0.1212121212121212;
        	const double inv = 0.03125;

#ifdef BBURST
        	double f,df;
        	double tcoast = TSTART;
        	if (t>=TSTART) {
        		const double bct = t0 + t1 + t2 + t3;
        		const double phase = std::floor((t-TSTART)/bct);
        		tcoast = TSTART + phase*bct;
        	}
        	double tfreeze = tcoast + t0;
        	double tburst = tfreeze + t1;
        	double tswim = tburst + t2;

        	if (t<tcoast) {
        		f = 1.0;
        		df = 0.0;
        	} else if (t<tfreeze) {
        		const double d = (t-tcoast)/(tfreeze-tcoast);
        		f = 1 - 3*d*d + 2*d*d*d;
        		df = 6*(d*d - d)/(tfreeze-tcoast);
        	} else if (t<tburst) {
        		f = 0.0;
        		df = 0.0;
        	} else if (t<tswim) {
        		const double d = (t-tburst)/(tswim-tburst);
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

        void _computeMidlineCoordinates(const double time)
        {
        	const double rampFac = rampFactorSine(time, Tperiod);
        	rX[0] = 0.0;
        	rY[0] =   rampFac*midline(rS[0], time, length, Tperiod, phaseShift);

        	for(int i=1;i<Nm;++i) {
        		rY[i]=rampFac*midline(rS[i], time, length, Tperiod, phaseShift);
        		const double dy = rY[i]-rY[i-1];
        		const double ds = rS[i] - rS[i-1];
        		const double dx = std::sqrt(ds*ds-dy*dy);
        		rX[i] = rX[i-1] + dx;
        	}
        }

        void _computeMidlineVelocities(const double time)
        {
        	const double rampFac =    rampFactorSine(time, Tperiod);
        	const double rampFacVel = rampFactorVelSine(time, Tperiod);

        	vX[0] = 0.0; //rX[0] is constant
        	vY[0] = rampFac*midlineVel(rS[0],time,length,Tperiod, phaseShift) +
        			rampFacVel*midline(rS[0],time,length,Tperiod, phaseShift);

        	for(int i=1;i<Nm;++i) {
        		vY[i]=rampFac*midlineVel(rS[i],time,length,Tperiod, phaseShift) +
        			  rampFacVel*midline(rS[i],time,length,Tperiod, phaseShift);
        		const double dy = rY[i]-rY[i-1];
        		const double dx = rX[i]-rX[i-1];
        		const double dVy = vY[i]-vY[i-1];
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

        void computeMidline(const double time)
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
    
    double volume_internal, J_internal, CoM_internal[2], vCoM_internal[2], theta_internal, angvel_internal, angvel_internal_prev, CoM_interpolated[3];
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

	void save(const int step_id, const double t, std::string filename = std::string()) override;
	void restart(const double t, std::string filename = std::string()) override;
    void update(const int step_id, const double t, const double dt, const double *Uinf) override;
    void getCenterOfMass(double CM[3]) const override;
    void create(const int step_id,const double time, const double dt, const double *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
