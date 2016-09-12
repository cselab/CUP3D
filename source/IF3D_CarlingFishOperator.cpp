//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_CarlingFishOperator.h"
#include <array>

namespace Fish
{
    inline double width(const double s, const Real L)
    {
    	if(s<0 or s>L) return 0;

    	const double sb=0.04*L;
    	const double st=0.95*L;
    	const double wt=0.01*L;
    	const double wh=0.04*L;

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
        const int Nm;
        double linMom[2], vol, J, angMom; // for diagnostics
        // start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
        const int iFishStart, iFishEnd;

    private:
        const Real time;
        const Real length;
        const Real Tperiod;
        const Real phaseShift;
        double Rmatrix2D[2][2];
        double Rmatrix3D[3][3];

        inline void _rotate2D(double &x, double &y)
        {
        	const double p[2] = {x,y};
        	x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
        	y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
        }

        inline void _translateAndRotate2D(const double pos[2], double &x, double &y)
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

        inline void _subtractAngularVelocity(const double angvel, const double x, const double y, double & vx, double & vy)
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
        
        inline double _integrationFac1(const int idx)
        {
            return width[idx]*height[idx];
        }
        
        inline double _integrationFac2(const int idx)
        {
            const double dnorXi = _d_ds(idx, norX, Nm);
            const double dnorYi = _d_ds(idx, norY, Nm);
            
            return 0.25*std::pow(width[idx],3)*height[idx]*(dnorXi*norY[idx] - dnorYi*norX[idx]);
        }
        
        inline double _integrationFac3(const int idx)
        {
            const double drXi = _d_ds(idx, rX, Nm);
            const double drYi = _d_ds(idx, rY, Nm);
           
//            return 0.25*std::pow(width[idx],3)*height[idx]*(drXi*norY[idx] - drYi*norX[idx]);
            return 0.25*std::pow(width[idx],3)*height[idx];
        }
        
        inline void _updateVolumeIntegration(const double fac1, const double ds, double & vol)
        {
            vol+=0.5*fac1*ds;
        }

        inline void _updateCoMIntegration(const int idx, const double fac1, const double fac2, const double ds, double CoM[2])
        {
            CoM[0] += 0.5*(rX[idx]*fac1 + norX[idx]*fac2)*ds;
            CoM[1] += 0.5*(rY[idx]*fac1 + norY[idx]*fac2)*ds;
        }
        
        inline void _updateLinMomIntegration(const int idx, const double fac1, const double fac2, const double ds, double linMom[2])
        {
            linMom[0] += 0.5*(vX[idx]*fac1 + vNorX[idx]*fac2)*ds;
            linMom[1] += 0.5*(vY[idx]*fac1 + vNorY[idx]*fac2)*ds;
        }
        
        inline void _updateAngMomIntegration(const int idx, const double fac1, const double fac2, const double fac3, const double ds, double & angMom)
        {
            double tmpSum = 0.0;
            tmpSum += (rX[idx]*vY[idx] - rY[idx]*vX[idx])*fac1;
            tmpSum += (rX[idx]*vNorY[idx] - rY[idx]*vNorX[idx] + vY[idx]*norX[idx] - vX[idx]*norY[idx])*fac2;
            tmpSum += (norX[idx]*vNorY[idx] - norY[idx]*vNorX[idx])*fac3;
            angMom += 0.5*tmpSum*ds;
        }
        
        inline void _updateJIntegration(const int idx, const double fac1, const double fac2, const double fac3, const double ds, double & J)
        {
            double tmpSum = 0.0;
            tmpSum += (rX[idx]*rX[idx] + rY[idx]*rY[idx])*fac1;
            tmpSum += 2.0*(rX[idx]*norX[idx] + rY[idx]*norY[idx])*fac2;
//            tmpSum += (norX[idx]*norX[idx]+norY[idx]*norY[idx])*fac3;
            tmpSum += fac3;
            J += 0.5*tmpSum*ds;
        }
        
public:
        
        FishMidlineData(const int Nm, const Real time, const Real length, const Real Tperiod, const Real phaseShift):
        	Nm(Nm),time(time),length(length),Tperiod(Tperiod),phaseShift(phaseShift),rS(_alloc(Nm)),rX(_alloc(Nm)),rY(_alloc(Nm)),
			vX(_alloc(Nm)),vY(_alloc(Nm)),norX(_alloc(Nm)),norY(_alloc(Nm)),vNorX(_alloc(Nm)),vNorY(_alloc(Nm)),
			width(_alloc(Nm)),height(_alloc(Nm)),iFishStart(0),iFishEnd(Nm-1)
		{
        	// fill s without extension
        	for(int i=0;i<Nm;++i)
        		rS[i] = length * 0.5 * (1.0 - std::cos(i * M_PI/((double)Nm-1))); // cosine: more points near head and tail
        	// rS[i] = i*length/((double)Nm-1); // linear: equally distributed points
		}

        FishMidlineData(const int N, const Real t, const Real len, const Real Tp, const Real phase, const std::pair<int,Real> ext_info):
        	Nm(N),time(t),length(len),Tperiod(Tp),phaseShift(phase),rS(_alloc(Nm)),rX(_alloc(Nm)),rY(_alloc(Nm)),
			vX(_alloc(Nm)),vY(_alloc(Nm)),norX(_alloc(Nm)),norY(_alloc(Nm)),vNorX(_alloc(Nm)),vNorY(_alloc(Nm)),
			width(_alloc(Nm)),height(_alloc(Nm)),iFishStart(ext_info.first),iFishEnd(Nm-1-ext_info.first)
        {
        	// extension_info contains number of extension points and extension dx
        	const int Next = extension_info.first; // number of points per extension
        	const int Nint = Nm - 2*Next; // number of interior points

        	// extension head
        	for(int i=0;i<Next;++i)
        		rS[i] = 0.0 - (Next- i) * extension_info.second;
        	// interior points
        	for(int i=0;i<Nint;++i)
        		rS[i+Next] = length * 0.5 * (1.0 - std::cos(i * M_PI/((double)Nint-1))); // cosine: more points near head and tail
        	// rS[i] = i*length/((double)Nint-1); // linear: equally distributed points
        	// extension tail
        	for(int i=0;i<Next;++i)
        		rS[i+Nint+Next] = length + (i + 1)*extension_info.second;
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
        
        void _correctTrajectory(const double dtheta, const double time, const double dt) {}

        void execute(const double time, const double l_tnext, const vector<Real>& input) {}

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

            for(int i=1;i<Nm-1;++i) {
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
            
            for(int i=0;i<Nm;++i) {
                _rotate2D(rX[i],rY[i]);
                _rotate2D(vX[i],vY[i]);
                _subtractAngularVelocity(angvel_internal, rX[i], rY[i], vX[i], vY[i]);           
            }

            _computeMidlineNormals();
        }
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
    	void _computeMidline()
    	{
    		const double _1oL = 1./lenght;
    		const std::array<double ,6> curvature_values = {
    				0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
					3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
    		};
    		const std::array<double ,6> curvature_points = {
    			0., .15*lenght, .4*lenght, .65*lenght, .9*lenght, lenght
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
    		const double _1oL = 1./length;
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
            _computeWidthsHeights();
    	}

    public:
    	CurvatureDefinedFishData(const int Nm, const Real t, const Real length, const Real Tperiod,
    			const Real phaseShift, const bool cosine=true)
    	: FishMidlineData(Nm,t,length,Tperiod,phaseShift,cosine), l_Tp(Tperiod), timeshift(0), time0(0),
		rK(_alloc(Nm)),vK(_alloc(Nm)),rC(_alloc(Nm)),vC(_alloc(Nm)),rA(_alloc(Nm)),vA(_alloc(Nm)),rB(_alloc(Nm)),vB(_alloc(Nm))
    	{
    	}

    	CurvatureDefinedFishData(const int Nm, const Real t, const Real length, const Real Tperiod,
    			const Real phaseShift, const std::pair<int,Real> extension_info, const bool cosine=true)
    	: FishMidlineData(Nm,t,length,Tperiod,phaseShift,extension_info,cosine), l_Tp(Tperiod), timeshift(0), time0(0),
		rK(_alloc(Nm)),vK(_alloc(Nm)),rC(_alloc(Nm)),vC(_alloc(Nm)),rA(_alloc(Nm)),vA(_alloc(Nm)),rB(_alloc(Nm)),vB(_alloc(Nm))
    	{
    	}

    	void _correctTrajectory(const double dtheta, const double time, const double dt) override
    	{
    		tVec tmp_curv;
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
    		_dealloc(rA);
    		_dealloc(vA);
    		_dealloc(rB);
    		_dealloc(vC);
    	}
    };

    struct CarlingFishMidlineData : public FishMidlineData
    {
    	inline double rampFactorSine(const Real t, const Real T)
    	{
    		return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
    	}

    	inline double rampFactorVelSine(const Real t, const Real T)
    	{
    		return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
    	}

    	inline double midline(const double s, const Real t, const Real L, const Real T, const Real phaseShift)
    	{
    		const double arg = 2.0*M_PI*(s/L - t/T + phaseShift);
    		return .4/3.3 *  (s + 0.03125*L)*std::sin(arg);
    	}

    	inline double midlineVel(const double s, const Real t, const Real L, const Real T, const Real phaseShift)
    	{
    		const double arg = 2.0*M_PI*(s/L - t/T + phaseShift);
    		return - .4/3.3 * (s + 0.03125*L) * (2.0*M_PI/T) * std::cos(arg);
    	}

        void _computeMidlineCoordinates()
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

        void _computeMidlineVelocities()
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

        void _computeMidlineNormals()
        {
            for(int i=0; i<Nm-1; i++) {
                const double ds = rS[i+1]-rS[i];
                const double tX = rX[i+1]-rX[i];
                const double tY = rY[i+1]-rY[i];
                norX[i] = -tY/ds;
                norY[i] =  tX/ds;

                const double tVX = vX[i+1]-vX[i];
                const double tVY = vY[i+1]-vY[i];

                vNorX[i] = -tVY/ds;
                vNorY[i] =  tVX/ds;
            }
            norX[Nm-1] = norX[Nm-2];
            norY[Nm-1] = norY[Nm-2];

            vNorX[Nm-1] = vNorX[Nm-2];
            vNorY[Nm-1] = vNorY[Nm-2];
        }

        void computeMidline()
        {
            _computeMidlineCoordinates();
            _computeMidlineVelocities();
            _computeMidlineNormals();
            _computeWidthsHeights();
        }
    };
}

struct VolumeSegment_OBB
{
    std::pair<int, int> s_range;
    double normalI[3]; // should be normalized and >=0
    double normalJ[3];
    double normalK[3];
    double w[3]; // halfwidth
    double c[3]; // center
    
    VolumeSegment_OBB(std::pair<int, int> s_range, const double bbox[3][2]):s_range(s_range)
    {
        normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
        normalI[0]=normalJ[1]=normalK[2]=1.0;
        for(int i=0;i<3;++i) {
            w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
            c[i] = bbox[i][0] + w[i];
            assert(w[i]>0);
        }
    }
    
    void normalizeNormals()
    {
        const double magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
        const double magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
        const double magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
        assert(magI > std::numeric_limits<double>::epsilon());
        assert(magJ > std::numeric_limits<double>::epsilon());
        assert(magK > std::numeric_limits<double>::epsilon());        
        const double invMagI = 1.0/magI;
        const double invMagJ = 1.0/magJ;
        const double invMagK = 1.0/magK;
        
        for(int i=0;i<3;++i) {
	    // also take absolute value since thats what we need when doing intersection checks later
            normalI[i]=std::abs(normalI[i])*invMagI;
            normalJ[i]=std::abs(normalJ[i])*invMagJ;
            normalK[i]=std::abs(normalK[i])*invMagK;
        }
    }

    void changeToComputationalFrame(const double position[3], const double quaternion[4])
    {
        // we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
        
        const double w = quaternion[0];
        const double x = quaternion[1];
        const double y = quaternion[2];
        const double z = quaternion[3];
        
        const double Rmatrix3D[3][3] = {
            {1-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
            {2*(x*y+z*w),    1-2*(x*x+z*z),  2*(y*z-x*w)},
            {2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x*x+y*y)}
        };
        
        const double p[3] = {c[0],c[1],c[2]};
                
        const double nx[3] = {normalI[0],normalI[1],normalI[2]};
        const double ny[3] = {normalJ[0],normalJ[1],normalJ[2]};
        const double nz[3] = {normalK[0],normalK[1],normalK[2]};
        
        for(int i=0;i<3;++i) {
            c[i] = Rmatrix3D[i][0]*p[0] + Rmatrix3D[i][1]*p[1] + Rmatrix3D[i][2]*p[2];
            
            normalI[i] = Rmatrix3D[i][0]*nx[0] + Rmatrix3D[i][1]*nx[1] + Rmatrix3D[i][2]*nx[2];
            normalJ[i] = Rmatrix3D[i][0]*ny[0] + Rmatrix3D[i][1]*ny[1] + Rmatrix3D[i][2]*ny[2];
            normalK[i] = Rmatrix3D[i][0]*nz[0] + Rmatrix3D[i][1]*nz[1] + Rmatrix3D[i][2]*nz[2];
        }
        
        c[0]+=position[0];
        c[1]+=position[1];
        c[2]+=position[2];
        
        normalizeNormals();
    }
    
    bool isIntersectingWithAABB(const double start[3],const double end[3], const double safe_distance = 0.0) const
    {
        const double AABB_w[3] = {
            0.5*(end[0] - start[0]) + 2.0*safe_distance,
            0.5*(end[1] - start[1]) + 2.0*safe_distance,
            0.5*(end[2] - start[2]) + 2.0*safe_distance
        }; // halfwidth
        
        const double AABB_c[3] = {
            start[0] + AABB_w[0] - safe_distance,
            start[1] + AABB_w[1] - safe_distance,
            start[2] + AABB_w[2] - safe_distance
        };
        
        assert(AABB_w[0]>0);
        assert(AABB_w[1]>0);
        assert(AABB_w[2]>0);
        
        bool intersects = true;
        double r;
        {
            r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
            intersects &= ((c[0]-r <= AABB_c[0] + AABB_w[0]) && (c[0]+r >= AABB_c[0] - AABB_w[0]));
            
            r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
            intersects &= ((c[1]-r <= AABB_c[1] + AABB_w[1]) && (c[1]+r >= AABB_c[1] - AABB_w[1]));
            
            r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
            intersects &= ((c[2]-r <= AABB_c[2] + AABB_w[2]) && (c[2]+r >= AABB_c[2] - AABB_w[2]));
        }
        {
            r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
            intersects &= ((AABB_c[0]-r <= c[0] + w[0]) && (AABB_c[0]+r >= c[0] - w[0]));
            
            r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
            intersects &= ((AABB_c[1]-r <= c[1] + w[1]) && (AABB_c[1]+r >= c[1] - w[1]));
            
            r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
            intersects &= ((AABB_c[2]-r <= c[2] + w[2]) && (AABB_c[2]+r >= c[2] - w[2]));
        }
        return intersects;
    }
};

struct PutFishOnBlocks
{    
    const CarlingFish::CarlingFishMidlineData & cfish;
    double position[3];
    double quaternion[4];
    double Rmatrix3D[3][3];

    PutFishOnBlocks(const CarlingFish::CarlingFishMidlineData & cfish, const double pos[3], const double quat[4]):
	cfish(cfish)
	{
        position[0]=pos[0];
        position[1]=pos[1];
        position[2]=pos[2];
        quaternion[0]=quat[0];
        quaternion[1]=quat[1];
        quaternion[2]=quat[2];
        quaternion[3]=quat[3];
        
        computeRotationMatrix();
	}

    PutFishOnBlocks(const PutFishOnBlocks& c):
	//segmentsPerBlock(c.segmentsPerBlock),deformation_velocities(c.deformation_velocities),
	cfish(c.cfish)
	{
        position[0]=c.position[0];
        position[1]=c.position[1];
        position[2]=c.position[2];
        quaternion[0]=c.quaternion[0];
        quaternion[1]=c.quaternion[1];
        quaternion[2]=c.quaternion[2];
        quaternion[3]=c.quaternion[3];
        
        computeRotationMatrix();        
	}
    
    void computeRotationMatrix()
    {
        const double w = quaternion[0];
        const double x = quaternion[1];
        const double y = quaternion[2];
        const double z = quaternion[3];
        
        const double R[3][3] = {
            {1-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
            {2*(x*y+z*w),    1-2*(x*x+z*z),  2*(y*z-x*w)},
            {2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x*x+y*y)}
        };
        
        memcpy(Rmatrix3D, R, sizeof(R));
    }
    
    inline int find_closest_dist(const int s, const int dir, const double x[3], double & oldDistSq) const
    {        
        if((s+dir)<cfish.iFishStart or (s+dir)>cfish.iFishEnd)
            return s;
        
        const double newDistSq = (x[0]-cfish.rX[s+dir])*(x[0]-cfish.rX[s+dir]) + (x[1]-cfish.rY[s+dir])*(x[1]-cfish.rY[s+dir]) + (x[2])*(x[2]);
        
        if(oldDistSq<=newDistSq)
            return s;
        else {
            oldDistSq = newDistSq;
            return s+dir;
        }
    }
    
    void changeVelocityToComputationalFrame(double x[3]) const
    {
        const double p[3] = {x[0],x[1],x[2]};
        
        // rotate (around CoM)
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
        x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
        x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];
    }
    
    template<typename T>
    void changeToComputationalFrame(T x[3]) const
    {
        const T p[3] = {x[0],x[1],x[2]};
        
        // rotate (around CoM)
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
        x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
        x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];
        
        // translate
        x[0]+=position[0];
        x[1]+=position[1];
        x[2]+=position[2];
    }
    
    template<typename T>
    void changeFromComputationalFrame(T x[3]) const
    {
        // translate back to CoM
        const T p[3] = {
            x[0]-position[0],
            x[1]-position[1],
            x[2]-position[2]
        };
        
        // rotate back around CoM
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[1][0]*p[1] + Rmatrix3D[2][0]*p[2];
        x[1]=Rmatrix3D[0][1]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[2][1]*p[2];
        x[2]=Rmatrix3D[0][2]*p[0] + Rmatrix3D[1][2]*p[1] + Rmatrix3D[2][2]*p[2];
    }

    double getSmallerDistToMidline(const int start_s, const double x[3], int & final_s) const
    {
        double relX[3] = {x[0],x[1],x[2]};
        changeFromComputationalFrame(relX);
        
        const double curDistSq = std::pow(relX[0]-cfish.rX[start_s],2) + std::pow(relX[1]-cfish.rY[start_s],2) + std::pow(relX[2],2);
        
        double distSq;
        // check right
        distSq = curDistSq;
        const int sRight = find_closest_dist(start_s, +1, relX, distSq);
        // check left
        distSq = curDistSq;
        const int sLeft = find_closest_dist(start_s, -1, relX, distSq);

        if(sRight==start_s and sLeft==start_s) {
            final_s = start_s;
            return distSq;
        }
        
        assert(sRight==start_s or sLeft==start_s);
        
        int curr_s = start_s;
        int new_s = sRight == start_s ? sLeft : sRight;
        const int dir = new_s-curr_s;
        while(curr_s not_eq new_s) {
            curr_s = new_s;
            new_s = find_closest_dist(curr_s,dir,relX, distSq);
        }

        final_s = new_s;
        return distSq;
    }

	//inline void operator()(const BlockInfo& info, FluidBlock3D& b) const
	void operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
	{
        Real org[3];
        info.pos(org, 0, 0, 0);        
        const Real invh = 1.0/info.h_gridpoint;
        
        {
            FluidElement * const e = &b(0,0,0);
            Real * const chi = &(defblock->chi[0][0][0]);
            Real * const udef = &(defblock->udef[0][0][0][0]);
            
            static const int n = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
            for(int i=0; i<n; i++) {
                e[i].tmp = 0; //reset psi here because were gonna do plus-equal next
                chi[i]=-1;
                udef[3*i+0]=0;
                udef[3*i+1]=0;
                udef[3*i+2]=0;
            }
        }
        
        // construct the shape (P2M with min(distance) as kernel)
        
        for(int i=0;i<vSegments.size();++i) {
        	for(int ss=std::max(vSegments[i].s_range.first,cfish.iFishStart); ss<=std::min(vSegments[i].s_range.second, cfish.iFishEnd); ++ss)
        	{
        		assert(ss>=cfish.iFishStart && ss<=cfish.iFishEnd);

        		// fill chi
        		const Real offset = cfish.height[ss] > cfish.width[ss] ? 0.5*M_PI : 0.0; // assume width is major axis, else correction

        		const Real ell_a = (Real)std::max(cfish.height[ss],cfish.width[ss]);
        		//                    const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : 0.25*info.h[0]/ell_a;
        		// maximum distance between two points is 2*ell_a * sin(theta). set this distance to dx/2 -->
        		const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : std::abs(std::asin(0.5*(0.25*info.h[0])/ell_a));

        		const int Ntheta = (int)std::ceil(2.0*M_PI/dtheta_target) + 1;
        		const Real dtheta = 2.0*M_PI/((Real)Ntheta);

        		for(int tt=0;tt<Ntheta;++tt) {
        			const Real theta = tt*dtheta + offset;
        			// create a surface point
        			// special treatment of tail (width = 0 --> no ellipse, just line)
        			const double hght = cfish.width[ss] == 0 ? cfish.height[ss]*(2*tt/((Real)Ntheta-1) - 1)
        					: cfish.height[ss]*std::sin(theta);
        			double myP[3] = {
        					(cfish.rX[ss] + cfish.width[ss]*std::cos(theta)*cfish.norX[ss]),
							(cfish.rY[ss] + cfish.width[ss]*std::cos(theta)*cfish.norY[ss]),
							hght
        			};

        			changeToComputationalFrame(myP);

        			const int iap[3] = {
        					(int)std::floor((myP[0]-org[0])*invh),
							(int)std::floor((myP[1]-org[1])*invh),
							(int)std::floor((myP[2]-org[2])*invh)
        			};

        			// support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
					const int start[3] = {
							std::max(-1, 0 - iap[0] ),
							std::max(-1, 0 - iap[1] ),
							std::max(-1, 0 - iap[2] )
					};
					const int end[3] = {
							std::min(+3, B::sizeX - iap[0]),
							std::min(+3, B::sizeY - iap[1]),
							std::min(+3, B::sizeZ - iap[2])
					};

					const double myP_distToMidlineSq = cfish.width[ss] == 0 ? std::pow(hght,2) :
							(double)(std::pow(cfish.width[ss]*std::cos(theta),2) + std::pow(cfish.height[ss]*std::sin(theta),2));


					if(myP_distToMidlineSq<std::numeric_limits<double>::epsilon()) {
						// if ss==iFishStart or ss==iFishEnd, the surface point lies on the midline --> myP_distToMidlineSq=0.
								// here our in/out criterion fails
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<B::sizeX);
							assert(idx[1]>=0 && idx[1]<B::sizeY);
							assert(idx[2]>=0 && idx[2]<B::sizeZ);

							double p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const double diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const double distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const double distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);

							changeFromComputationalFrame(p);
							const double distPlanar = std::sqrt( std::pow(p[0]-cfish.rX[closest_s],2) +
									std::pow(p[1]-cfish.rY[closest_s],2) );
							const double distHeight = std::abs(p[2]);
							const Real sign = (distPlanar > cfish.width[closest_s] or distHeight > cfish.height[closest_s]) ? -1.0 : 1.0;

							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					} else {
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<B::sizeX);
							assert(idx[1]>=0 && idx[1]<B::sizeY);
							assert(idx[2]>=0 && idx[2]<B::sizeZ);

							double p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const double diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const double distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const double distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);

							const Real sign = distToMidlineSq >= myP_distToMidlineSq ? -1.0 : 1.0;

							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					}
        		}
        	}
        }
        
        // construct the deformation velocities (P2M with hat function as kernel)
        
        for(int i=0;i<vSegments.size();++i) {
        	for(int ss=vSegments[i].s_range.first;ss<=vSegments[i].s_range.second;++ss) {
        		assert(ss>=0 && ss<=cfish.Nm-1);
        		// P2M udef of a slice at this s
        		const double myWidth = (ss < cfish.iFishStart ? cfish.width[cfish.iFishStart]
													: (ss > cfish.iFishEnd ? cfish.width[cfish.iFishEnd] : cfish.width[ss]));
        		const double myHeight = (ss < cfish.iFishStart ? cfish.height[cfish.iFishStart]
													: (ss > cfish.iFishEnd ? cfish.height[cfish.iFishEnd] : cfish.height[ss]));
        		const double ds_defGrid = info.h_gridpoint;
        		// towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice (and so we can floor!)
        		const double extension = 3*info.h_gridpoint;

        		const int Nh = std::floor( (myHeight+extension)/ds_defGrid );

        		for(int ih=-Nh;ih<=Nh; ++ih) {
        			const Real offsetH = ih*ds_defGrid;
        			// add an extra extension when width == 0 (to deal with large curvatures near head and/or tail):
        			const double currentWidth = myWidth== 0 ? extension : myWidth * std::sqrt(1 - std::pow(offsetH/(myHeight+extension),2));
        			const double actualWidth = (cfish.height[ss] == 0 or std::abs(offsetH)>=cfish.height[ss]) ? 0.0
        										: cfish.width[ss] * std::sqrt(1 - std::pow(offsetH/cfish.height[ss],2));
        			const int Nw = std::floor( (currentWidth+extension)/ds_defGrid); // add extension here to make sure we have it in each direction

        			for(int iw=-Nw;iw<=Nw; ++iw) {
        				const Real offsetW = iw*ds_defGrid;

        				double xp[3] = {
        						(cfish.rX[ss] + offsetW*cfish.norX[ss]),
								(cfish.rY[ss] + offsetW*cfish.norY[ss]),
								offsetH
        				};

        				changeToComputationalFrame(xp);

        				xp[0] = (xp[0]-org[0])*invh;
        				xp[1] = (xp[1]-org[1])*invh;
        				xp[2] = (xp[2]-org[2])*invh;

        				double udef[3] = {
        						(cfish.vX[ss] + offsetW*cfish.vNorX[ss]),
								(cfish.vY[ss] + offsetW*cfish.vNorY[ss]),
								0.0
        				};

        				changeVelocityToComputationalFrame(udef);

        				const Real ap[3] = {
        						std::floor((Real)xp[0]),
								std::floor((Real)xp[1]),
								std::floor((Real)xp[2])
        				};

        				const int iap[3] = {
        						(int)ap[0],
								(int)ap[1],
								(int)ap[2]
        				};

        				// now we P2M
						const int start[3] = {
								std::max(0, 0 - iap[0] ),
								std::max(0, 0 - iap[1] ),
								std::max(0, 0 - iap[2] )
        				};
        				const int end[3] = {
        						std::min(+2, B::sizeX - iap[0]),
								std::min(+2, B::sizeY - iap[1]),
								std::min(+2, B::sizeZ - iap[2])
        				};

        				Real wghts[3][2];
        				for(int c=0;c<3;++c) {
        					const Real t[2] = {
        							std::abs((Real)xp[c] - (ap[c]+0)),
									std::abs((Real)xp[c] - (ap[c]+1))
        					};
        					wghts[c][0] = 1.0 - t[0];
        					wghts[c][1] = 1.0 - t[1];
        				}

        				const bool isInside = (std::abs(offsetW) < actualWidth) && (std::abs(offsetH) < cfish.height[ss]);

        				for(int sz=start[2]; sz<end[2];++sz) {
        					const Real wz = wghts[2][sz];
        					for(int sy=start[1];sy<end[1];++sy) {
        						const Real wywz = wz*wghts[1][sy];
        						for(int sx=start[0];sx<end[0];++sx) {
        							const Real wxwywz = wywz*wghts[0][sx];
        							assert(wxwywz>=0 && wxwywz<=1);

        							const int idx[3] = {
        									iap[0] + sx,
											iap[1] + sy,
											iap[2] + sz,
        							};
        							assert(idx[0]>=0 && idx[0]<B::sizeX);
        							assert(idx[1]>=0 && idx[1]<B::sizeY);
        							assert(idx[2]>=0 && idx[2]<B::sizeZ);

        							defblock->udef[idx[2]][idx[1]][idx[0]][0] += wxwywz*udef[0];
        							defblock->udef[idx[2]][idx[1]][idx[0]][1] += wxwywz*udef[1];
        							defblock->udef[idx[2]][idx[1]][idx[0]][2] += wxwywz*udef[2];

        							b(idx[0],idx[1],idx[2]).tmp += wxwywz;

        							// set sign for all interior points
									if( (std::abs(defblock->chi[idx[2]][idx[1]][idx[0]] + 1) < 5*std::numeric_limits<Real>::epsilon()) && isInside)
										defblock->chi[idx[2]][idx[1]][idx[0]] = 1.0;
        						}
        					}
        				}
        			}
        		}
        	}
        }
        
        // finalize signed distance function in psi
        {
            for(int iz=0; iz<FluidBlock3D::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock3D::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock3D::sizeX; ix++) {
                        const Real normfac = b(ix,iy,iz).tmp > std::numeric_limits<Real>::epsilon() ? b(ix,iy,iz).tmp : 1.0;
                        defblock->udef[iz][iy][ix][0] /= normfac;
                        defblock->udef[iz][iy][ix][1] /= normfac;
                        defblock->udef[iz][iy][ix][2] /= normfac;
			// change from signed squared distance function to normal sdf
                        b(ix,iy,iz).tmp = (defblock->chi[iz][iy][ix] > 0 ?
                        		std::sqrt(defblock->chi[iz][iy][ix]) : -std::sqrt(-defblock->chi[iz][iy][ix]));
                        
                    }
        }
	}
};

struct PutFishOnBlocks_Finalize
{
    Real t;
	int stencil_start[3], stencil_end[3];

    PutFishOnBlocks_Finalize(const double dummy): t(0)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}
    
    PutFishOnBlocks_Finalize(const PutFishOnBlocks_Finalize& c): t(c.t)
	{
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
    }
    
	template<typename Lab>
	void operator()(Lab & lab,const BlockInfo& info,FluidBlock& b,ObstacleBlock* const defblock,surfaceBlocks& surf,double& tmp[4]) const
	{
		const double h = info.h_gridpoint;
		const double inv2h = 0.5/h;
		const double invh2 = 1/(h*h);
		const double eps = std::numeric_limits<Real>::epsilon();

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++)
		{
			double p[3];
			info.pos(p, ix,iy,iz);
			/*
				const bool in_band =
				lab(ix,iy,iz)*lab(ix+1,iy,iz) <=0 ||
				lab(ix,iy,iz)*lab(ix-1,iy,iz) <=0 ||
				lab(ix,iy,iz)*lab(ix,iy+1,iz) <=0 ||
				lab(ix,iy,iz)*lab(ix,iy-1,iz) <=0 ||
				lab(ix,iy,iz)*lab(ix,iy,iz+1) <=0 ||
				lab(ix,iy,iz)*lab(ix,iy,iz-1) <=0;
			 */
			if (lab(ix,iy).tmp >= +3*h || lab(ix,iy).tmp <= -3*h) {
				//if(not in_band or (std::abs(defblock->chi[iz][iy][ix])>0.5)) {
				const double H = lab(ix,iy,iz).tmp > 0 ? 1.0 : 0.0;
				tmp[0]+=H;
				tmp[1]+=p[0]*H;
				tmp[2]+=p[1]*H;
				tmp[3]+=p[2]*H;

				defblock->chi[iz][iy][ix] = H;
				b(ix,iy,iz).tmp = std::max(defblock->chi[iz][iy][ix], b(ix,iy,iz).tmp);
				continue;
			}

			const double distPx = lab(ix+1,iy,iz).tmp;
			const double distMx = lab(ix-1,iy,iz).tmp;
			const double distPy = lab(ix,iy+1,iz).tmp;
			const double distMy = lab(ix,iy-1,iz).tmp;
			const double distPz = lab(ix,iy,iz+1).tmp;
			const double distMz = lab(ix,iy,iz-1).tmp;

			const double IplusX = distPx < 0 ? 0 : distPx;
			const double IminuX = distMx < 0 ? 0 : distMx;
			const double IplusY = distPy < 0 ? 0 : distPy;
			const double IminuY = distMy < 0 ? 0 : distMy;
			const double IplusZ = distPz < 0 ? 0 : distPz;
			const double IminuZ = distMz < 0 ? 0 : distMz;

			const double HplusX = distPx == 0 ? 0.5 : (distPx < 0 ? 0 : 1);
			const double HminuX = distMx == 0 ? 0.5 : (distMx < 0 ? 0 : 1);
			const double HplusY = distPy == 0 ? 0.5 : (distPy < 0 ? 0 : 1);
			const double HminuY = distMy == 0 ? 0.5 : (distMy < 0 ? 0 : 1);
			const double HplusZ = distPz == 0 ? 0.5 : (distPz < 0 ? 0 : 1);
			const double HminuZ = distMz == 0 ? 0.5 : (distMz < 0 ? 0 : 1);

			// gradU
			const double gradUX = inv2h * (distPx - distMx);
			const double gradUY = inv2h * (distPy - distMy);
			const double gradUZ = inv2h * (distPz - distMz);
			const double gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ;

			// gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
			const double gradIX = inv2h * (IplusX - IminuX);
			const double gradIY = inv2h * (IplusY - IminuY);
			const double gradIZ = inv2h * (IplusZ - IminuZ);
			const double numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;

			const double gradHX = inv2h * (HplusX - HminuX);
			const double gradHY = inv2h * (HplusY - HminuY);
			const double gradHZ = inv2h * (HplusZ - HminuZ);
			const double numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;

			const double Delta = std::abs(gradUSq) < eps ? numD : numD/gradUSq;
			const double H     = std::abs(gradUSq) < eps ? numH : numH/gradUSq;

			if (Delta>1e-6) {
				const double dchidx = -Delta*gradUX;
				const double dchidy = -Delta*gradUY;
				const double dchidz = -Delta*gradUZ;
				surf.add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, Delta);
			}
			tmp[0]+=H;
			tmp[1]+=p[0]*H;
			tmp[2]+=p[1]*H;
			tmp[3]+=p[2]*H;

			defblock->chi[iz][iy][ix] = H;
			b(ix,iy,iz).tmp = std::max(defblock->chi[iz][iy][ix], b(ix,iy,iz).tmp);
		}
	}
};

void IF3D_CarlingFishOperator::create()
{
    // STRATEGY
    // we need some things already
    // - the internal angle at the previous timestep, obtained from integrating the actual deformation velocities
	// 						 (not the imposed deformation velocies, because they dont have zero angular momentum)
    // - the internal angular velocity at previous timestep
    
    // 1. create midline
    // 2. integrate to find CoM, angular velocity, etc
    // 3. shift midline to CoM frame: zero internal linear momentum and angular momentum
    
    // 4. split the fish into segments (according to s)
    // 5. rotate the segments to computational frame (comp CoM and angle)
    // 6. for each Block in the domain, find those segments that intersect it
    // 7. for each of those blocks, allocate an ObstacleBlock
    
    // 8. put the 3D shape on the grid: SDF-P2M for sdf, normal P2M for udef
    // 9. create the Chi out of the SDF. In same sweep, compute the actual CoM
    // 10. compute all shit: linear momentum, angular momentum etc.
    // 11. correct deformation velocity to nullify momenta for the final discrete representation
    
    const Real min_dx = (1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel());
    const Real target_ds = 0.5*min_dx;
    const Real target_Nm = length/target_ds;
  
    // multiple of 100: TODO why?
    const int Nm = 100*(int)std::ceil(target_Nm/100) + 1;
//    const int Nm = 1201;
    const int Nsegments = 100;
    assert((Nm-1)%Nsegments==0);

    // deal with extension
    const Real dx_extension = 0.25*min_dx;
    const int Nextension = 12;// up to 3dx on each side (to get proper interpolation up to 2dx)
    // 1.
    myFish.computeMidline();
    
    // 2. & 3.
    volume_internal = myFish.integrateLinearMomentum(CoM_internal, vCoM_internal);
    assert(volume_internal > std::numeric_limits<double>::epsilon());
    myFish.changeToCoMFrameLinear(CoM_internal, vCoM_internal);
    
    angvel_internal_prev = angvel_internal;
    myFish.integrateAngularMomentum(angvel_internal);
    J_internal = myFish.J;
    // update theta now with new angvel info
    //theta_internal -= 0.5*sim_dt*(angvel_internal + angvel_internal_prev); // negative: we subtracted this angvel
    myFish.changeToCoMFrameAngular(theta_internal, angvel_internal);

#ifndef DNDEBUG
    {
        double dummy_CoM_internal[2], dummy_vCoM_internal[2], dummy_angvel_internal;

        // check that things are zero
        const double volume_internal_check = myFish.integrateLinearMomentum(dummy_CoM_internal,dummy_vCoM_internal);
        myFish.integrateAngularMomentum(dummy_angvel_internal);
        
        assert(std::abs(dummy_CoM_internal[0])<std::numeric_limits<Real>::epsilon());
        assert(std::abs(dummy_CoM_internal[1])<std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish.linMom[0])<std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish.linMom[1])<std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish.angMom)<std::numeric_limits<Real>::epsilon());
        assert(std::abs(volume_internal - volume_internal_check) < std::numeric_limits<Real>::epsilon());
    }
#endif
    
    // 4.
    std::vector<VolumeSegment_OBB> vSegments;
    {
        for(int i=0;i<Nsegments;++i)
        {
            const int idx = i * (Nm-1)/Nsegments;
            const int next_idx = (i+1) * (Nm-1)/Nsegments;
            
            // find bounding box based on this
            double bbox[3][2]={
                {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()},
                {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()},
                {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()}
            };
            for(int ss=idx; ss<=next_idx; ++ss) {
                const double xBnd[2] = {(myFish.rX[ss] - myFish.norX[ss]*myFish.width[ss]), (myFish.rX[ss] + myFish.norX[ss]*myFish.width[ss])};
                const double yBnd[2] = {(myFish.rY[ss] - myFish.norY[ss]*myFish.width[ss]), (myFish.rY[ss] + myFish.norY[ss]*myFish.width[ss])};
                const double zBnd[2] = {(-myFish.height[ss]), (+myFish.height[ss])};
                bbox[0][0] = std::min({bbox[0][0],xBnd[0],xBnd[1]});
                bbox[0][1] = std::max({bbox[0][1],xBnd[0],xBnd[1]});
                bbox[1][0] = std::min({bbox[1][0],yBnd[0],yBnd[1]});
                bbox[1][1] = std::max({bbox[1][1],yBnd[0],yBnd[1]});
                bbox[2][0] = std::min({bbox[2][0],zBnd[0],zBnd[1]});
                bbox[2][1] = std::max({bbox[2][1],zBnd[0],zBnd[1]});
            }
            
            // create a new segment
            VolumeSegment_OBB volumeSegment(std::make_pair(idx, next_idx), bbox);
            vSegments.push_back(volumeSegment);
        }
    }
    assert(vSegments.size()==Nsegments);

    // 5.
    for(int i=0;i<Nsegments;++i)
        vSegments[i].changeToComputationalFrame(position,quaternion);
    
    // clear deformation velocities
    for(auto & entry : obstacleBlocks)
        delete entry.second;
    obstacleBlocks.clear();

    // 6. & 7.    
    std::map<int, std::vector<VolumeSegment_OBB> > segmentsPerBlock;
    {
        vector<BlockInfo> vInfo = grid.getBlocksInfo();
        for(int i=0;i<vInfo.size();++i) {
            const BlockInfo & info = vInfo[i];
            double pStart[3], pEnd[3];
            info.pos(pStart, 0, 0, 0);
            info.pos(pEnd, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            const double safe_distance = 2.0*info.h_gridpoint; // two points on each side

            for(int s=0;s<Nsegments;++s)
                if(vSegments[s].isIntersectingWithAABB(pStart,pEnd,safe_distance))
                    segmentsPerBlock[i].push_back(vSegments[s]);
            
            // allocate new blocks if necessary
            if(segmentsPerBlock.find(i) != segmentsPerBlock.end()) {
                assert(obstacleBlocks.find(i) == obstacleBlocks.end());
                obstacleBlocks[i] = new ObstacleBlock;
                obstacleBlocks[i]->clear();
            }
        }
    }
    
    assert(not segmentsPerBlock.empty());
    assert(segmentsPerBlock.size() == deformation_velocities.size());
    
    // 8.
    {
#pragma omp parallel
    	{
        	PutFishOnBlocks putfish(myFish, position, quaternion);

#pragma omp for schedule(static)
    		for(int i=0; i<vInfo.size(); i++) {
    			auto pos = segmentsPerBlock.find(i);
    			FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;

    			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					b(ix,iy,iz).tmp = -1; // negative unity

    			if(pos == segmentsPerBlock.end()) continue;
    			assert(obstacleBlocks.find(i) != obstacleBlocks.end());
    			ObstacleBlock* const defblock = obstacleBlocks.find(i)->second;

    			putfish(ary[i], b, defblock, pos->second);
    		}
    	}
    }
    
    // 9. & 10. & 11.
    {
    	const int nthreads = omp_get_max_threads();
    	vector<surfaceBlocks> dataPerThread(nthreads);

    	double sumVX(0), sumVY(0), sumVZ(0), sumXs(0);
#pragma omp parallel
    	{
    		Lab lab;
    		lab.prepare(*grid, stencil_start,  stencil_end, true);
    		FishObstacle::PutFishOnBlocks_Finalize finalize(0.);
    		const int tid = omp_get_thread_num();

#pragma omp for schedule(static), reduction(+:sumVX,sumVY,sumVZ,sumXs)
    		for(int i=0; i<vInfo.size(); i++) {
    			auto pos = segmentsPerBlock.find(i);
    			if(pos == segmentsPerBlock.end()) continue;
    			FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;
    			assert(obstacleBlocks.find(i) != obstacleBlocks.end());
    			ObstacleBlock* const defblock = obstacleBlocks.find(i)->second;

    			lab.load(ary[i], 0);
    			double tmp[4];
    			finalize(lab, ary[i], b, defblock, pos->second, dataPerThread(tid), tmp);
    			sumVX+=tmp[0];
    			sumVY+=tmp[1];
    			sumVZ+=tmp[2];
    			sumXs+=tmp[3];
    		}
    	}

    	double totVX(0), totVY(0), totVZ(0), totXs(0);
    	MPI::COMM_WORLD.Allreduce(&sumVX, &totVX, 1, MPI::DOUBLE, MPI::SUM);
    	MPI::COMM_WORLD.Allreduce(&sumVY, &totVY, 1, MPI::DOUBLE, MPI::SUM);
    	MPI::COMM_WORLD.Allreduce(&sumVZ, &totVZ, 1, MPI::DOUBLE, MPI::SUM);
    	MPI::COMM_WORLD.Allreduce(&sumXs, &totXs, 1, MPI::DOUBLE, MPI::SUM);

    	surfData.finalizeOnGrid(dataPerThread);

        assert(totXs > std::numeric_limits<double>::epsilon());
        CoM_interpolated[0]=totVX/totXs;
        CoM_interpolated[1]=totVY/totXs;
        CoM_interpolated[2]=totVZ/totXs;

        _makeDefVelocitiesMomentumFree(CoM_interpolated);

        /*
#ifndef DNDEBUG
        {
            ComputeAll computeall(vInfo,grid.getBlockCollection(),CoM_interpolated,deformation_velocities);
            tbb::parallel_reduce(blocked_range<int>(0,vInfo.size(),1),computeall);
            
            std::cout << computeall.properties.linearMomentum[0] << std::endl;
            std::cout << computeall.properties.linearMomentum[1] << std::endl;
            std::cout << computeall.properties.linearMomentum[2] << std::endl;
            std::cout << computeall.properties.angularMomentum[0] << std::endl;
            std::cout << computeall.properties.angularMomentum[1] << std::endl;
            std::cout << computeall.properties.angularMomentum[2] << std::endl;
            
            assert(std::abs(computeall.properties.linearMomentum[0]) < std::numeric_limits<Real>::epsilon());
            assert(std::abs(computeall.properties.linearMomentum[1]) < std::numeric_limits<Real>::epsilon());
            assert(std::abs(computeall.properties.linearMomentum[2]) < std::numeric_limits<Real>::epsilon());
            assert(std::abs(computeall.properties.angularMomentum[0]) < std::numeric_limits<Real>::epsilon());
            assert(std::abs(computeall.properties.angularMomentum[1]) < std::numeric_limits<Real>::epsilon());
            assert(std::abs(computeall.properties.angularMomentum[2]) < std::numeric_limits<Real>::epsilon());
        }
#endif
        */
    }
    
    //bObstacleBlocksFilled=true;
}

void IF3D_CarlingFishOperator::update(const double t, const double dt)
{
    // synchronize internal time
    sim_time = t + dt;
    sim_dt = dt;
    // update position and angles
    IF3D_ObstacleOperator::update(t, dt);
    // negative: we subtracted this angvel
    theta_internal -= sim_dt*angvel_internal;
    // tell solver we have to redo the fish
    //bObstacleBlocksFilled = false;
}

void IF3D_CarlingFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF2D_ObstacleOperator::_parseArguments(parser);

	parser.set_strict_mode();
	length = parser("-length").asDouble();
	Tperiod = parser("-T").asDouble();
	phaseShift = parser("-phi").asDouble();

	bCorrectTrajectory = parser("-Correct").asBool(false);
    Tstartlearn = parser("-Tstartlearn").asDouble(1e6);
    GoalDX = parser("-GoalDX").asDouble(0.0);
    nActions = parser("-nActions").asInt(0);
	/*
    randomStart = parser("-randomStart").asBool(false);
    if (randomStart) {
    	printf("Random start\n");
    	std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<Real> dis(-1.,1.);
    	const Real rseed1 = .8*length*dis(gen);
    	const Real rseed2 = .2*length*dis(gen);
    	const Real rseed3 = .1* M_PI *dis(gen);
    	position[0] += rseed1/sqrt(2.)-rseed2/sqrt(2.);
    	position[1] += rseed1/sqrt(2.)+rseed2/sqrt(2.);
    	ext_pos[0] = position[0];
    	ext_pos[1] = position[1];
    	angle += rseed3;
    }
    sr->updateInstant(position[0], ext_pos[0], position[1], ext_pos[1], angle, 0., 0., 0.);
    //TODO state and reward:
    sr->t_next_comm = Tstartlearn - 1/2.; //i want to reset time-averages before first actual comm
    bool bForgiving = parser("-easyFailBox").asBool(false);
    sr->bForgiving = bForgiving;
    sr->GoalDX = GoalDX;
    sr->thExp = angle;

    randomActions = parser("-randomActions").asBool(false);
    if (randomActions) printf("Fish doing random turns\n");
    useLoadedActions = parser("-useLoadedActions").asBool(false);
    if (useLoadedActions) {
        Real dummy_time;
        vector<Real> action(nActions);
        ifstream in("orders_1.txt"); //FUCKING TODO NEED TO USE SOME POINTERS IN THIS SHIT
        std::string line;
        if(in.good()) {
            while (getline(in, line)) {
                istringstream line_in(line);
                line_in >> dummy_time;
                line_in >> action[0];
                if(nActions==2) line_in >> action[1];
                //i want to do pop back later:
                loadedActions.insert(loadedActions.begin(),action);
            }
        } else { printf("Could not load actions from file orders_1.txt\n"); abort(); }
        in.close();
    }
    */
}

/*
void IF2D_StefanLearnTurnOperator::execute(Communicator * comm, const int iAgent, const double time)
{
    if (time < Tstartlearn) {
        sr->resetAverage();
        sr->t_next_comm = Tstartlearn;

        //TMP: first rnd action is going to be taken after a while
        if (randomActions) sr->t_next_comm = Tstartlearn+6.;

        return;
    }

    if (not bInteractive) {
        if (not randomActions) { sr->t_next_comm=1e6; return; }
        //we might decide to pause turning, then just pick a pause counter
        if (nPauseActions-- > 0) {
            vector<Real> raction(1,0.);
            myFish->execute(time, sr->t_next_comm, raction);
            sr->t_next_comm += .5*myFish->Tperiod;
            printf("pausing from turning at time %f, will still pause %d turns. Next turn at time %f\n",time,nPauseActions,sr->t_next_comm);
            return;
        }
        vector<Real> raction(1,signLastTurn);
        myFish->execute(time, sr->t_next_comm, raction);
        sr->t_next_comm += .5*myFish->Tperiod;
        printf("turning at time %f with modifier %f, turn counter is %d. Next turn at time %f\n",time,signLastTurn,nTurnActions,sr->t_next_comm);

        if (++nTurnActions >= 4) {
            nPauseActions = 10.;
            nTurnActions = 0;
            signLastTurn *= -1.;
        }

    } else if (useLoadedActions) {

        vector<Real> actions(nActions);
        if (loadedActions.size()>1) {
            actions = loadedActions.back();
            loadedActions.pop_back();
        } //else zero actions
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }
        sr->resetAverage();

    } else {

        const Real relT= fmod(time,1.); //1 is Tperiod
#ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*NpLatLine : 25+  8*NpLatLine;
#else
        const int nStates = (nActions==1) ? 20+10*NpLatLine : 25+ 10*NpLatLine;
#endif
        vector<Real> state(nStates), actions(nActions);

        int k(0);
        state[k++] = sr->Xrel - GoalDX;
        state[k++] = sr->Yrel;
        state[k++] = sr->RelAng;
        state[k++] = relT;
        state[k++] = new_curv;
        state[k++] = old_curv;

        if(nActions==2) { //this is for backwards compatibility
            state[k++] = new_Tp;
                        //2.*M_PI*((time-time0)/l_Tp +timeshift -rS[i]/length) + M_PI*phaseShift
            Real Fshift = 2.*((-myFish->time0)/myFish->l_Tp +myFish->timeshift)+myFish->phaseShift;
            Fshift = fmod(Fshift,2.0);
            state[k++] = (Fshift<0) ? 2.+Fshift : Fshift;
            state[k++] = sr->VX;
            state[k++] = sr->VY;
            state[k++] = sr->AV;
        }

        state[k++] = sr->Dist;
        state[k++] = sr->Quad;
        state[k++] = sr->VxAvg;
        state[k++] = sr->VyAvg;
        state[k++] = sr->AvAvg;
        state[k++] = sr->Pout;
        state[k++] = sr->defPower;
        state[k++] = sr->EffPDef;
        state[k++] = sr->PoutBnd;
        state[k++] = sr->defPowerBnd;
        state[k++] = sr->EffPDefBnd;
        state[k++] = sr->Pthrust;
        state[k++] = sr->Pdrag;
        state[k++] = sr->ToD;

        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVBelow[j];
        #ifndef _NOVISION_
        for (int j=0; j<2*NpLatLine; j++) state[k++] = sr->raySight[j];
        #endif
        const Real reward = (sr->info==2) ? -10 : sr->EffPDefBnd;
        comm->sendState(iAgent-1, sr->info, state, reward); //TODO
        fflush(0);
        if (sr->info==2) return;

        sr->info = 0;

        comm->recvAction(actions);
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }

        #ifndef TRAINING
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        if(nActions==2)
            filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
        #endif //TRAINING

        sr->resetAverage();
    }
}
*/

void IF3D_CarlingFishOperator::getCenterOfMass(double& CM[3]) const
{
	// return computation CoM, not the one were advecting
	CM[0]=CoM_interpolated[0];
	CM[1]=CoM_interpolated[1];
	CM[2]=CoM_interpolated[2];
}

void IF3D_CarlingFishOperator::save(const double t, string filename)
{
    assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
    
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<double>::digits10 + 1);
    
    if(filename==std::string())
        savestream.open("restart_IF3D_Carling.txt");
    else
        savestream.open(filename+".txt");

    savestream << sim_time << "\t" << sim_dt << std::endl;
    savestream << position[0] << "\t" << position[1] << "\t" << position[2] << std::endl;
    savestream << quaternion[0] << "\t" << quaternion[1] << "\t" << quaternion[2] << "\t" << quaternion[3] << std::endl;
    savestream << transVel[0] << "\t" << transVel[1] << "\t" << transVel[2] << std::endl;
    savestream << angVel[0] << "\t" << angVel[1] << "\t" << angVel[2] << std::endl;
    savestream << theta_internal << "\t" << angvel_internal << std::endl;    
    savestream.close();
    
}

void IF3D_CarlingFishOperator::restart(const double t, string filename)
{
    std::ifstream restartstream;
    
    if(filename==std::string())
        restartstream.open("restart_IF3D_Carling.txt");
    else
        restartstream.open(filename+".txt");
    
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal;
    restartstream.close();
    
    {
        std::cout << "RESTARTED FISH: " << std::endl;
        std::cout << "TIME, DT: \t" << sim_time << " " << sim_dt << std::endl;
        std::cout << "POS: \t" << position[0] << " " << position[1] << " " << position[2] << std::endl;
        std::cout << "ANGLE: \t" << quaternion[0] << "\t" << quaternion[1] << "\t" << quaternion[2] << "\t" << quaternion[3] << std::endl;
        std::cout << "TVEL: \t" << transVel[0] << " " << transVel[1] << " " << transVel[2] << std::endl;
        std::cout << "AVEL: \t" << angVel[0] << " " << angVel[1] << " " << angVel[2] << std::endl;
        std::cout << "INTERN: \t" << theta_internal << " " << angvel_internal << std::endl;
    }
}

void IF3D_CarlingFishOperator::writeToFile(const int step_id, const Real t, std::string filename)
{
    std::string fname = (filename==std::string()) ? "fish" : filename;

    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<double>::digits10 + 1);
    
    savestream.open(fname+"_interpolated.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";        
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "CoM[2]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "linMom[2]\t";
        savestream << "angMom[0]\t";
        savestream << "angMom[1]\t";
        savestream << "angMom[2]" << std::endl;
    }    
        
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << object_ongrid.volume << "\t";
    savestream << CoM_interpolated[0] << "\t";
    savestream << CoM_interpolated[1] << "\t";
    savestream << CoM_interpolated[2] << "\t";
    savestream << object_ongrid.linearMomentum[0] << "\t";
    savestream << object_ongrid.linearMomentum[1] << "\t";
    savestream << object_ongrid.linearMomentum[2] << "\t";
    savestream << object_ongrid.angularMomentum[0] << "\t";
    savestream << object_ongrid.angularMomentum[1] << "\t";
    savestream << object_ongrid.angularMomentum[2] << "\t";
    savestream << object_ongrid.J[2] << std::endl;
    savestream.close();
    
    savestream.open(fname+"_internal.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "angMom\t";
        savestream << "theta\t";
        savestream << "angvel" << std::endl;
    }
    
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << volume_internal << "\t";
    savestream << CoM_internal[0] << "\t";
    savestream << CoM_internal[1] << "\t";
    savestream << vCoM_internal[0]*volume_internal << "\t";
    savestream << vCoM_internal[1]*volume_internal << "\t";
    savestream << angvel_internal*J_internal << "\t";
    savestream << theta_internal << "\t";
    savestream << angvel_internal << std::endl;
    savestream.close();
    
    savestream.open(fname+"_computation.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "pos[0]\t";
        savestream << "pos[1]\t";
        savestream << "pos[2]\t";
        savestream << "quat[0]\t";
        savestream << "quat[1]\t";
        savestream << "quat[2]\t";
        savestream << "quat[3]\t";
        savestream << "transVel[0]\t";
        savestream << "transVel[1]\t";
        savestream << "transVel[2]\t";
        savestream << "angVel[0]\t";
        savestream << "angVel[1]\t";
        savestream << "angVel[2]" << std::endl;        
    }
    
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << position[0] << "\t";
    savestream << position[1] << "\t";
    savestream << position[2] << "\t";
    savestream << quaternion[0] << "\t";
    savestream << quaternion[1] << "\t";
    savestream << quaternion[2] << "\t";
    savestream << quaternion[3] << "\t";
    savestream << transVel[0] << "\t";
    savestream << transVel[1] << "\t";
    savestream << transVel[2] << "\t";
    savestream << angVel[0] << "\t";
    savestream << angVel[1] << "\t";
    savestream << angVel[2] << std:: endl;
    savestream.close();
}
