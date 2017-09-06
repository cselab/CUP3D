//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_FishLibrary__
#define __IncompressibleFluids3D__IF3D_FishLibrary__
//#define BBURST
#include <cmath>
#include <array>
#include "IF2D_Interpolation1D.h"
#include "GenericOperator.h"
#include "IF2D_Frenet.h"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>

const int NPPEXT = 3; //was 3
const int TGTPPB = 2; //was 2 i think
const int NEXTDX = 2;
#define __BSPLINE

namespace Schedulers
{
template<int Npoints>
struct ParameterScheduler
{
	std::array<Real, Npoints> parameters_t0; // parameters at t0
	std::array<Real, Npoints> parameters_t1; // parameters at t1
	std::array<Real, Npoints> dparameters_t0; // derivative at t0
	Real t0, t1; // t0 and t1

	void save(std::string filename)
	{
		std::ofstream savestream;
		savestream.setf(std::ios::scientific);
		savestream.precision(std::numeric_limits<Real>::digits10 + 1);
		savestream.open(filename+".txt");

		savestream << t0 << "\t" << t1 << std::endl;
		for(int i=0;i<Npoints;++i)
			savestream << parameters_t0[i]  << "\t"
								 << parameters_t1[i]  << "\t"
								 << dparameters_t0[i] << std::endl;
		savestream.close();
	}

	void restart(std::string filename)
	{
		std::ifstream restartstream;
		restartstream.open(filename+".txt");
		int rank;
                MPI_Comm_rank(MPI_COMM_WORLD,&rank);
		if(!rank) std::cout << filename << " ";
		restartstream >> t0 >> t1;
		for(int i=0;i<Npoints;++i) {
			restartstream >> parameters_t0[i] >> parameters_t1[i] >> dparameters_t0[i];
			if(!rank) 
			std::cout<<parameters_t0[i]<<" "<<parameters_t1[i]<<" "<<dparameters_t0[i];
		}
		if(!rank) std::cout << std::endl;
		restartstream.close();
	}

	ParameterScheduler()
	{
		t0=-1; t1=0;
		parameters_t0 = std::array<Real, Npoints>();
		parameters_t1 = std::array<Real, Npoints>();
		dparameters_t0 = std::array<Real, Npoints>();
	}

	void transition(const Real t, const Real tstart, const Real tend,
			const std::array<Real, Npoints> parameters_tend,
			const bool UseCurrentDerivative = false)
	{
		if(t<tstart or t>tend) return; // this transition is out of scope
		//if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

		// we transition from whatever state we are in to a new state
		// the start point is where we are now: lets find out
		std::array<Real, Npoints> parameters;
		std::array<Real, Npoints> dparameters;
		gimmeValues(tstart,parameters,dparameters);

		/*
            if (Npoints == 7)
                printf("[t0 t1] were [%f %f], will be [%f %f], parameters %f dpdt %f param end %f\n",
                		t0,t1,tstart,tend,parameters[3],dparameters[3], parameters_tend[3]);
		 */

		// fill my members
		t0 = tstart;
		t1 = tend;
		parameters_t0 = parameters;
		parameters_t1 = parameters_tend;
		dparameters_t0 = UseCurrentDerivative ? dparameters : std::array<Real, Npoints>();
	}

	void transition(const Real t, const Real tstart, const Real tend,
			const std::array<Real, Npoints> parameters_tstart,
			const std::array<Real, Npoints> parameters_tend)
	{
		if(t<tstart or t>tend) return; // this transition is out of scope
		if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

		// fill my members
		t0 = tstart;
		t1 = tend;
		parameters_t0 = parameters_tstart;
		parameters_t1 = parameters_tend;
	}

	void gimmeValues(const Real t,std::array<Real, Npoints> & parameters, std::array<Real, Npoints> & dparameters)
	{
		// look at the different cases
		if(t<t0 or t0<0) { // no transition, we are in state 0
			parameters = parameters_t0;
			dparameters = std::array<Real, Npoints>();
		} else if(t>t1) { // no transition, we are in state 1
			parameters = parameters_t1;
			dparameters = std::array<Real, Npoints>();
		} else { // we are within transition: interpolate
			for(int i=0;i<Npoints;++i)
				IF2D_Interpolation1D::cubicInterpolation(t0,t1,t,parameters_t0[i],parameters_t1[i],dparameters_t0[i],0.0,parameters[i],dparameters[i]);
		}
	}

	void gimmeValues(const Real t,std::array<Real, Npoints> & parameters)
	{
		std::array<Real, Npoints> dparameters_whocares; // no derivative info
		return gimmeValues(t,parameters,dparameters_whocares);
	}
};

struct ParameterSchedulerScalar : ParameterScheduler<1>
{
	void transition(const Real t, const Real tstart, const Real tend, const Real parameter_tend, const bool UseCurrentDerivative = false)
	{
		const std::array<Real, 1> myParameter = {parameter_tend};
		return ParameterScheduler<1>::transition(t,tstart,tend,myParameter,UseCurrentDerivative);
	}

	void gimmeValues(const Real t, Real & parameter, Real & dparameter)
	{
		std::array<Real, 1> myParameter, mydParameter;
		ParameterScheduler<1>::gimmeValues(t, myParameter, mydParameter);
		parameter = myParameter[0];
		dparameter = mydParameter[0];
	}

	void gimmeValues(const Real t, Real & parameter)
	{
		std::array<Real, 1> myParameter;
		ParameterScheduler<1>::gimmeValues(t, myParameter);
		parameter = myParameter[0];
	}
};

template<int Npoints>
struct ParameterSchedulerVector : ParameterScheduler<Npoints>
{
	void gimmeValues(const Real t, const std::array<Real, Npoints> & positions, const int Nfine,
			const Real * const positions_fine, Real * const parameters_fine, Real * const dparameters_fine)
	{
		// we interpolate in space the start and end point
		Real parameters_t0_fine[Nfine];
		Real parameters_t1_fine[Nfine];
		Real dparameters_t0_fine[Nfine];

		IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->parameters_t0.data(), Npoints, positions_fine, parameters_t0_fine,  Nfine);
		IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->parameters_t1.data(), Npoints, positions_fine, parameters_t1_fine,  Nfine);
		IF2D_Interpolation1D::naturalCubicSpline(positions.data(), this->dparameters_t0.data(),Npoints, positions_fine, dparameters_t0_fine, Nfine);

		// look at the different cases
		if(t<this->t0 or this->t0<0) {
			// no transition, we are in state 0
			for(int i=0;i<Nfine;++i) {
				parameters_fine[i] = parameters_t0_fine[i];
				dparameters_fine[i] = 0.0;
			}
		} else if(t>this->t1) {
			// no transition, we are in state 1
			for(int i=0;i<Nfine;++i) {
				parameters_fine[i] = parameters_t1_fine[i];
				dparameters_fine[i] = 0.0;
			}
		} else {
			// we are within transition: interpolate in time for each point of the fine discretization
			for(int i=0;i<Nfine;++i)
				IF2D_Interpolation1D::cubicInterpolation(this->t0,this->t1,t,parameters_t0_fine[i],parameters_t1_fine[i],dparameters_t0_fine[i],
						0.0,				   parameters_fine[i],	 dparameters_fine[i]);
		}
	}

	void gimmeValues(const Real t,std::array<Real, Npoints> & parameters)
	{
		ParameterScheduler<Npoints>::gimmeValues(t, parameters);
	}

	void gimmeValues(const Real t,std::array<Real, Npoints> & parameters, std::array<Real, Npoints> & dparameters)
	{
		ParameterScheduler<Npoints>::gimmeValues(t, parameters, dparameters);
	}
};

template<int Npoints>
struct ParameterSchedulerLearnWave : ParameterScheduler<Npoints>
{
	void gimmeValues(const Real t, const Real Twave, const Real Length,
		const std::array<Real, Npoints> & positions, const int Nfine,
		const Real* const positions_fine, Real* const parameters_fine, Real* const dparameters_fine)
	{
		const Real _1oL = 1./Length;
		const Real _1oT = 1./Twave;
		// the fish goes through (as function of t and s) a wave function that describes the curvature
		for(int i=0;i<Nfine;++i) {
			const Real c = positions_fine[i]*_1oL - (t - this->t0)*_1oT; //traveling wave coord
			bool bCheck = true;

			if (c < positions[0]) { // Are you before latest wave node?
				IF2D_Interpolation1D::cubicInterpolation(
					c, positions[0], c,
					this->parameters_t0[0], this->parameters_t0[0],
					parameters_fine[i], dparameters_fine[i]);
				bCheck = false;
			}
			else if (c > positions[Npoints-1]) {// Are you after oldest wave node?
				IF2D_Interpolation1D::cubicInterpolation(
					positions[Npoints-1], c, c,
					this->parameters_t0[Npoints-1], this->parameters_t0[Npoints-1],
					parameters_fine[i], dparameters_fine[i]);
				bCheck = false;
			} else {
				for (int j=1; j<Npoints; ++j) { // Check at which point of the travelling wave we are
					if (( c >= positions[j-1] ) && ( c <= positions[j] )) {
						IF2D_Interpolation1D::cubicInterpolation(
							positions[j-1], positions[j], c,
							this->parameters_t0[j-1], this->parameters_t0[j],
							parameters_fine[i], dparameters_fine[i]);
						dparameters_fine[i] = -dparameters_fine[i]*_1oT; // df/dc * dc/dt
						bCheck = false;
					}
				}
			}
			if (bCheck) { std::cout << "Ciaone2!" << std::endl; abort(); }
		}
	}

	void Turn(const Real b, const Real t_turn) // each decision adds a node at the beginning of the wave (left, right, straight) and pops last node
	{
		this->t0 = t_turn;
		for(int i=Npoints-1; i>1; --i)
				this->parameters_t0[i] = this->parameters_t0[i-2];
		this->parameters_t0[1] = b;
		this->parameters_t0[0] = 0;
	}
};
}

struct FishSkin
{
    public:
    const int Npoints;
    Real * const xSurf;
    Real * const ySurf;
    Real * const normXSurf;
    Real * const normYSurf;
    Real * const midX;
    Real * const midY;

    FishSkin(const int N): Npoints(N), xSurf(_alloc(N)), ySurf(_alloc(N)),
    normXSurf(_alloc(N-1)), normYSurf(_alloc(N-1)), midX(_alloc(N-1)), midY(_alloc(N-1))
    {}

    virtual ~FishSkin()
    {
        _dealloc(xSurf);
        _dealloc(ySurf);
        _dealloc(normXSurf);
        _dealloc(normYSurf);
        _dealloc(midX);
        _dealloc(midY);
    }

    protected:
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
};

class FishMidlineData
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
	Real * const rXold;
	Real * const rYold;
	Real oldTime = 0.0;
	bool firstStep = true;
	const bool bDoubleHinge=false;
	Real sHinge2=0.0, torqueZsecMarkers=0.0, tOld=0.0, dTheta2=0.0;

	Real linMom[2], vol, J, angMom; // for diagnostics
	// start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
	const int iFishStart, iFishEnd;
	const Real length, Tperiod, phaseShift;
	Real l_Tp, time0, timeshift;
	Schedulers::ParameterSchedulerVector<6> curvScheduler;
	Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
	Schedulers::ParameterSchedulerVector<6> adjustScheduler;
	FishSkin * upperSkin, * lowerSkin;
	Real finSize = 1.1e-1, waveLength = 1.0;

 protected:
	Real Rmatrix2D[2][2];
	Real Rmatrix3D[3][3];

	#ifndef __BSPLINE
	inline Real _width(const Real s, const Real L)
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

	inline Real _height(const Real s, const Real L)
	{
		if(s<0 or s>L) return 0;
		const Real a=0.51*L;
		const Real b=0.08*L;
		return b*std::sqrt(1 - std::pow((s-a)/a,2));
	}
	#endif

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
			return 0.5*((vals[idx+1]-vals[idx])/(rS[idx+1]-rS[idx]) +
			            (vals[idx]-vals[idx-1])/(rS[idx]-rS[idx-1]) );
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

	void _prepareRotation2D(const Real angle)
	{
		Rmatrix2D[0][0] = Rmatrix2D[1][1] = std::cos(angle);
		Rmatrix2D[0][1] = -std::sin(angle);
		Rmatrix2D[1][0] = -Rmatrix2D[0][1];
	}

	#ifdef __BSPLINE
	/*
	function inputs: xc, yc are n sized arrays which contain the control points of the cubic b spline
	function outputs onto res: assumed to be either the width or the height
	*/
	void integrateBSpline(Real* const res, const Real* const xc,
																				 const Real* const yc, const int n)
	{
		double len = 0;
	  for (int i=0; i<n-1; i++) {
	    len += std::sqrt(std::pow(xc[i]-xc[i+1],2) +
	                     std::pow(yc[i]-yc[i+1],2));
	  }

		gsl_bspline_workspace *bw;
	  gsl_vector *B;
	  // allocate a cubic bspline workspace (k = 4)
	  bw = gsl_bspline_alloc(4, n-2);
	  B = gsl_vector_alloc(n);
	  gsl_bspline_knots_uniform(0.0, len, bw);

	  double ti = 0;
		for(int i=0;i<Nm;++i) {
			res[i] = 0;
			if (rS[i]>0 and rS[i]<length) {
				const double dtt = 0.1*(rS[i]-rS[i-1]);
				while (true) {
					double xi = 0;
			    gsl_bspline_eval(ti, B, bw);
			    for (int j=0; j<n; j++)
			      xi += xc[j]*gsl_vector_get(B, j);
					if (xi >= rS[i]) break;
					ti += dtt;
				}

				for (int j=0; j<n; j++)
					res[i] += yc[j]*gsl_vector_get(B, j);
			}
		}
  	gsl_bspline_free(bw);
  	gsl_vector_free(B);
	}
	#endif

	void _computeWidthsHeights()
	{
#ifdef __BSPLINE
		const int nh = 8;
		const Real xh[8] = {0, 0, .2*length, .4*length,
			.6*length, .8*length, length, length};
		// Slim Zebrafish
		//const Real yh[8] = {0, 5.5e-2*length, 6.8e-2*length, 7.6e-2*length,
		//	6.4e-2*length, 7.2e-3*length, 1.1e-1*length, 0};

		// Large fin
		//const Real yh[8] = {0, 5.5e-2*length, 1.8e-1*length, 2e-1*length,
		//	6.4e-2*length, 2e-3*length, 3.25e-1*length, 0};

		printf("TailFinSize = %f, Wavelength = %f\n", finSize, waveLength);
		fflush(NULL);
		const Real yh[8] = {0, 5.5e-2*length, 1.8e-1*length, 2e-1*length,
			6.4e-2*length, 2e-3*length, finSize*length, 0};

        const int nw = 6;
        const Real xw[6] = {0, 0, length/3., 2*length/3., length, length};
        const Real yw[6] = {0, 8.9e-2*length, 7.0e-2*length,
		3.0e-2*length, 2.0e-2*length, 0};
        //const Real yw[6] = {0, 8.9e-2*length, 1.7e-2*length,
	//	1.6e-2*length, 1.3e-2*length, 0};
		integrateBSpline(width,  xw, yw, nw);
		integrateBSpline(height, xh, yh, nh);
	  #else
		for(int i=0;i<Nm;++i) {
			width[i]  = _width(rS[i],length);
			height[i] = _height(rS[i],length);
		}
	  #endif
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

 public:
	FishMidlineData(const int Nm, const Real len, const Real Tp, const Real phase, const Real dx_ext):
		Nm(Nm),length(len),Tperiod(Tp),phaseShift(phase),l_Tp(Tperiod),timeshift(0),
		time0(0),rS(_alloc(Nm)),rX(_alloc(Nm)),rY(_alloc(Nm)),vX(_alloc(Nm)),
		vY(_alloc(Nm)),norX(_alloc(Nm)),norY(_alloc(Nm)),vNorX(_alloc(Nm)),
		vNorY(_alloc(Nm)),width(_alloc(Nm)),height(_alloc(Nm)),
		iFishStart(NEXTDX*NPPEXT),iFishEnd(Nm-1-NEXTDX*NPPEXT),
		rXold(_alloc(Nm)), rYold(_alloc(Nm))
	{

		std::fill(rXold, rXold+Nm, 0.0);
		std::fill(rYold, rYold+Nm, 0.0);

		// extension_info contains number of extension points and extension dx
		const int Nextension = NEXTDX*NPPEXT; // up to 3dx on each side (to get proper interpolation up to 2dx)
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
		// FinSize will not be read in from text file by Carling constructor. Call this again when need to over-write with updated values

    upperSkin = new FishSkin(Nint);
    lowerSkin = new FishSkin(Nint);
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
		_dealloc(rXold);
		_dealloc(rYold);
		if(upperSkin not_eq nullptr) {
			delete upperSkin;
			upperSkin=nullptr;
		}
		if(lowerSkin not_eq nullptr) {
			delete lowerSkin;
			lowerSkin=nullptr;
		}
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
			vX[i] += angvel_internal*rY[i];
			vY[i] -= angvel_internal*rX[i];
		}
		_computeMidlineNormals();
	}

	void computeSurface()
	{
	    const int Nskin = lowerSkin->Npoints;
	    // Compute surface points by adding width to the midline points
			#pragma omp parallel for
	    for(int i=0; i<Nskin; ++i)
	    {
	        Real normal[2] = {norX[iFishStart + i], norY[iFishStart + i]};
	        Real const norm_mod1 = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
	        normal[0] /= norm_mod1;
	        normal[1] /= norm_mod1;
					assert(width[iFishStart + i] >= 0);
	        lowerSkin->xSurf[i] = rX[iFishStart + i] - width[iFishStart + i]*normal[0];
	        lowerSkin->ySurf[i] = rY[iFishStart + i] - width[iFishStart + i]*normal[1];
	        upperSkin->xSurf[i] = rX[iFishStart + i] + width[iFishStart + i]*normal[0];
	        upperSkin->ySurf[i] = rY[iFishStart + i] + width[iFishStart + i]*normal[1];
	    }
	}

	void computeSkinNormals(const Real theta_comp, const Real CoM_comp[3])
	{
			_prepareRotation2D(theta_comp);
			for(int i=0; i<Nm; ++i) {
				_rotate2D(rX[i], rY[i]);
				rX[i] += CoM_comp[0];
				rY[i] += CoM_comp[1];
			}

			const int Nskin = lowerSkin->Npoints;
	    // Compute midpoints as they will be pressure targets
	    #pragma omp parallel for
	    for(int i=0; i<Nskin-1; ++i)
	    {
	        lowerSkin->midX[i] = (lowerSkin->xSurf[i] + lowerSkin->xSurf[i+1])/2.;
	        upperSkin->midX[i] = (upperSkin->xSurf[i] + upperSkin->xSurf[i+1])/2.;
	        lowerSkin->midY[i] = (lowerSkin->ySurf[i] + lowerSkin->ySurf[i+1])/2.;
	        upperSkin->midY[i] = (upperSkin->ySurf[i] + upperSkin->ySurf[i+1])/2.;

	        lowerSkin->normXSurf[i]=  (lowerSkin->ySurf[i+1]-lowerSkin->ySurf[i]);
	        upperSkin->normXSurf[i]=  (upperSkin->ySurf[i+1]-upperSkin->ySurf[i]);
	        lowerSkin->normYSurf[i]= -(lowerSkin->xSurf[i+1]-lowerSkin->xSurf[i]);
	        upperSkin->normYSurf[i]= -(upperSkin->xSurf[i+1]-upperSkin->xSurf[i]);

	        const Real normL = std::sqrt( std::pow(lowerSkin->normXSurf[i],2) +
																				std::pow(lowerSkin->normYSurf[i],2) );
					const Real normU = std::sqrt( std::pow(upperSkin->normXSurf[i],2) +
																				std::pow(upperSkin->normYSurf[i],2) );

	        lowerSkin->normXSurf[i] /= normL;
	        upperSkin->normXSurf[i] /= normU;
	        lowerSkin->normYSurf[i] /= normL;
	        upperSkin->normYSurf[i] /= normU;

	        //if too close to the head or tail, consider a point further in, so that we are pointing out for sure
	        const int ii = (i<8) ? 8 : ((i > Nskin-9) ? Nskin-9 : i);

	        const Real dirL = lowerSkin->normXSurf[i]*(lowerSkin->midX[i]-rX[iFishStart+ii])
	        								+ lowerSkin->normYSurf[i]*(lowerSkin->midY[i]-rY[iFishStart+ii]);
		      const Real dirU = upperSkin->normXSurf[i]*(upperSkin->midX[i]-rX[iFishStart+ii])
		      								+ upperSkin->normYSurf[i]*(upperSkin->midY[i]-rY[iFishStart+ii]);

	        if(dirL < 0) {
	            lowerSkin->normXSurf[i] *= -1.0;
	            lowerSkin->normYSurf[i] *= -1.0;
	        }
	        if(dirU < 0) {
	            upperSkin->normXSurf[i] *= -1.0;
	            upperSkin->normYSurf[i] *= -1.0;
	        }
	    }
	}

	void surfaceToCOMFrame(const double theta_internal, const double CoM_internal[2])
	{
	    _prepareRotation2D(theta_internal);
	    // Surface points rotation and translation

	    for(int i=0; i<upperSkin->Npoints; ++i)
	    //for(int i=0; i<upperSkin->Npoints-1; ++i)
	    {
	        upperSkin->xSurf[i] -= CoM_internal[0];
	        upperSkin->ySurf[i] -= CoM_internal[1];
	        _rotate2D(upperSkin->xSurf[i], upperSkin->ySurf[i]);
	        lowerSkin->xSurf[i] -= CoM_internal[0];
	        lowerSkin->ySurf[i] -= CoM_internal[1];
	        _rotate2D(lowerSkin->xSurf[i], lowerSkin->ySurf[i]);
	        #if 0
	        upperSkin->midX[i] -= CoM_internal[0];
	        upperSkin->midY[i] -= CoM_internal[1];
	        _rotate2D(upperSkin->midX[i], upperSkin->midY[i]);
	        lowerSkin->midX[i] -= CoM_internal[0];
	        lowerSkin->midY[i] -= CoM_internal[1];
	        _rotate2D(lowerSkin->midX[i], lowerSkin->midY[i]);

	        _rotate2D(upperSkin->normXSurf[i],upperSkin->normYSurf[i]);
	        _rotate2D(lowerSkin->normXSurf[i],lowerSkin->normYSurf[i]);
					#endif
	    }
			#if 0
	    upperSkin->xSurf[upperSkin->Npoints-1] -= CoM_internal[0];
	    upperSkin->ySurf[upperSkin->Npoints-1] -= CoM_internal[1];
	    _rotate2D(upperSkin->xSurf[upperSkin->Npoints-1],
	              upperSkin->ySurf[upperSkin->Npoints-1]);
	    lowerSkin->xSurf[lowerSkin->Npoints-1] -= CoM_internal[0];
	    lowerSkin->ySurf[lowerSkin->Npoints-1] -= CoM_internal[1];
	    _rotate2D(lowerSkin->xSurf[lowerSkin->Npoints-1],
	              lowerSkin->ySurf[lowerSkin->Npoints-1]);
			#endif
	}

	void surfaceToComputationalFrame(const double theta_comp, const double CoM_interpolated[3])
	{
	    _prepareRotation2D(theta_comp);

	    for(int i=0; i<upperSkin->Npoints; ++i)
	    //for(int i=0; i<upperSkin->Npoints-1; ++i)
	    {
	        _rotate2D(upperSkin->xSurf[i], upperSkin->ySurf[i]);
	        upperSkin->xSurf[i] += CoM_interpolated[0];
	        upperSkin->ySurf[i] += CoM_interpolated[1];
	        _rotate2D(lowerSkin->xSurf[i], lowerSkin->ySurf[i]);
	        lowerSkin->xSurf[i] += CoM_interpolated[0];
	        lowerSkin->ySurf[i] += CoM_interpolated[1];
	        #if 0
	        _rotate2D(upperSkin->midX[i], upperSkin->midY[i]);
	        upperSkin->midX[i] += CoM_interpolated[0];
	        upperSkin->midY[i] += CoM_interpolated[1];
	        _rotate2D(lowerSkin->midX[i], lowerSkin->midY[i]);
	        lowerSkin->midX[i] += CoM_interpolated[0];
	        lowerSkin->midY[i] += CoM_interpolated[1];

	        _rotate2D(upperSkin->normXSurf[i],upperSkin->normYSurf[i]);
	        _rotate2D(lowerSkin->normXSurf[i],lowerSkin->normYSurf[i]);
					#endif
	    }
			#if 0
	    _rotate2D(upperSkin->xSurf[upperSkin->Npoints-1],
	              upperSkin->ySurf[upperSkin->Npoints-1]);
	    upperSkin->xSurf[upperSkin->Npoints-1] += CoM_interpolated[0];
	    upperSkin->ySurf[upperSkin->Npoints-1] += CoM_interpolated[1];
	    _rotate2D(lowerSkin->xSurf[lowerSkin->Npoints-1],
	              lowerSkin->ySurf[lowerSkin->Npoints-1]);
	    lowerSkin->xSurf[lowerSkin->Npoints-1] += CoM_interpolated[0];
	    lowerSkin->ySurf[lowerSkin->Npoints-1] += CoM_interpolated[1];
			#endif
	}

	virtual void computeMidline(const Real time) = 0;
	virtual void _correctAmplitude(const Real dAmp, const Real vAmp, const Real time, const Real dt) {}
	virtual void _correctTrajectory(const Real dtheta, const Real vtheta, const Real time, const Real dt) {}
	virtual void execute(const Real time, const Real l_tnext, const vector<Real>& input) {}
};

class NacaMidlineData : public FishMidlineData
{
 protected:
	inline Real _naca_width(const double s, const Real L)
	{
			if(s<0 or s>L) return 0;
			const Real a = 0.2969;
			const Real b =-0.1260;
			const Real c =-0.3516;
			const Real d = 0.2843;
			const Real e =-0.1015;
			const Real t = 0.12*L;
			const Real p = s/L;
			return 5*t* (a*sqrt(p) +b*p +c*p*p +d*p*p*p + e*p*p*p*p);
	}

	void _naca_integrateBSpline(Real* const res, const Real* const xc,
													 const Real* const yc, const int n)
	{
		double len = 0;
	  for (int i=0; i<n-1; i++)
	    len += std::sqrt(std::pow(xc[i]-xc[i+1],2) +
	                     std::pow(yc[i]-yc[i+1],2));

		gsl_bspline_workspace *bw;
	  gsl_vector *B;
	  // allocate a cubic bspline workspace (k = 4)
	  bw = gsl_bspline_alloc(4, n-2);
	  B = gsl_vector_alloc(n);
	  gsl_bspline_knots_uniform(0.0, len, bw);

	  double ti = 0;
		for(int i=0;i<Nm;++i) {
			res[i] = 0;
			if (rS[i]>0 and rS[i]<length) {
				const double dtt = 0.1*(rS[i]-rS[i-1]);
				while (true) {
					double xi = 0;
			    gsl_bspline_eval(ti, B, bw);
			    for (int j=0; j<n; j++)
			      xi += xc[j]*gsl_vector_get(B, j);
					if (xi >= rS[i]) break;
					ti += dtt;
				}

				for (int j=0; j<n; j++)
					res[i] += yc[j]*gsl_vector_get(B, j);
			}
		}
  	gsl_bspline_free(bw);
  	gsl_vector_free(B);
	}

	void _naca_computeWidthsHeights()
	{
		/*
				const int nw = 7;
			  const Real xw[nw] = {0., 0., length*.25, length*.5, length*.75, length, length};
			  const Real yw[nw] = {0, .5*length, .5*length, .5*length, .5*length, .5*length, 0};
		*/
		const int nw = 5;
	  	const Real xw[nw] = {0., 0., length*.5, length, length};
	  	const Real yw[nw] = {0, .5*length, .5*length, .5*length, 0};

		_naca_integrateBSpline(width, xw, yw, nw);
		for(int i=0;i<Nm;++i)
			height[i]  = _naca_width(rS[i],length);
	}

 public:
	NacaMidlineData(const int Nm, const Real length, const Real dx_ext) :
	FishMidlineData(Nm,length,1,0,dx_ext)
	{
		//just overwrite default width and height
		_naca_computeWidthsHeights();
		rX[0] = rY[0] = vX[0] = vY[0] = 0;
		for(int i=1; i<Nm; ++i) {
			rY[i] = vX[i] = vY[i] = 0;
			rX[i] = rX[i-1] + std::fabs(rS[i]-rS[i-1]);
		}
		_computeMidlineNormals();
		#if 0
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);
			if (rank!=0) return;
			FILE * f = fopen("fish_profile","w");
	    for (int i=0; i<Nm; ++i) printf("%g %g %g %g %g\n",
			rX[i],rY[i],rS[i],width[i],height[i]);
			fclose(f);
			printf("Dumped midline\n");
		#endif
	}
	void computeMidline(const Real time) override
	{
		rX[0] = rY[0] = vX[0] = vY[0] = 0;
		for(int i=1; i<Nm; ++i) {
			rY[i] = vX[i] = vY[i] = 0;
			rX[i] = rX[i-1] + std::fabs(rS[i]-rS[i-1]);
		}
		_computeMidlineNormals();
	}
};


class CarlingFishMidlineData : public FishMidlineData
{
 protected:
	//burst-coast:
	const Real tStart;
	Real t0, t1, t2, t3, lowestAmp=1e12;
	const Real fac, inv;
	//const Real sHinge, AhingeTheta, ThingeTheta, hingePhi;
	const Real sHinge, ThingeTheta;
	Real AhingeTheta, hingePhi;
	const bool bBurst, bHinge;
	Real kSpring=0.0;
	const Real kMaxSpring=100.0; // Expect torque values on the order of 1e-5 at steady, and 1e-3 at startup
Real thetaOld = 0.0, avgTorque = 0.0, runningTorque = 0.0, timeNminus = 0.0;
int prevTransition = 0;
	const bool quadraticAmplitude;

	inline Real rampFactorSine(const Real t, const Real T) const
	{
		return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
	}

	inline Real rampFactorVelSine(const Real t, const Real T) const
	{
		return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
	}

	inline Real midline(const Real s, const Real t, const Real L, const Real T,
		const Real phaseShift) const
	{
		const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
		
		double yCurrent;
		if(quadraticAmplitude){
			yCurrent = (s*s*0.15/L) *std::sin(arg);
		} else {
			yCurrent = fac * (s + inv*L)*std::sin(arg);
		}

		if(bHinge){
			if(s>sHinge){
				double yNot;
				if(quadraticAmplitude){
					yNot =  (sHinge*sHinge*0.15/L)*std::sin(2.0*M_PI*(sHinge/(waveLength*L) - t/T + phaseShift));
				}else{
					yNot =  fac *  (sHinge + inv*L)*std::sin(2.0*M_PI*(sHinge/(waveLength*L) - t/T + phaseShift));
				}
				const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
				const double dydsNot = std::sin(currentTheta);
				yCurrent = yNot + dydsNot*(s-sHinge);
			}
                }
                return yCurrent;
	}

	inline Real midlineVel(const Real s, const Real t, const Real L, const Real T,
		const Real phaseShift) const
	{
		const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
		double velCurrent;
		if(quadraticAmplitude){
			velCurrent = - (s*s*0.15/L)*(2.0*M_PI/T)*std::cos(arg);
		}else{
			velCurrent = - fac*(s + inv*L)*(2.0*M_PI/T)*std::cos(arg);
		}

		if(bHinge){
			if(s>sHinge){
				//const double yNot =  4./33 *  (sHinge + 0.03125*L)*std::sin(2.0*M_PI*(sHinge/L - t/T + phaseShift));
				double velNot;
				if(quadraticAmplitude){
					velNot =  -2.0*M_PI/T * (sHinge*sHinge*0.15/L)*std::cos(2.0*M_PI*(sHinge/(L*waveLength) - t/T + phaseShift));
				}else{
					velNot =  -2.0*M_PI/T * fac *  (sHinge + inv*L)*std::cos(2.0*M_PI*(sHinge/(L*waveLength) - t/T + phaseShift));
				}
				const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
				const double currentThetaDot = AhingeTheta * 2.0*M_PI/ThingeTheta * std::cos(2.0*M_PI*(t/ThingeTheta + hingePhi));
				const double dydsNotDT = std::cos(currentTheta)*currentThetaDot;
				velCurrent = velNot + dydsNotDT*(s-sHinge);
			}
		}
		return velCurrent;
	}

	inline Real midlineBeC(const Real s, const Real t, const Real L, const Real T,
		const Real phaseShift, const Real f) const
	{
		const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
		return f * fac * (s + inv*L)*std::sin(arg);
	}

	inline Real midlineVelBeC(const Real s, const Real t, const Real L, const Real T,
		const Real phaseShift, const Real f, const Real df) const
	{
		const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
		return fac*(s + inv*L)*(df*std::sin(arg) - f*(2.0*M_PI/T)*std::cos(arg));
	}

        std::pair<double, double> cubicHermite(const double f1, const double f2, const double x){
                const double a = 2*(f1-f2);
                const double b = -3*(f1-f2);
                const double retVal = a*x*x*x + b*x*x + f1;
                const double deriv = 3*a*x*x + 2*b*x;
                return std::make_pair(retVal, deriv);
        }

	void _computeMidlineCoordinatesBeC(const Real time)
	{
		const Real rampFac = rampFactorSine(time, Tperiod);
		const Real rampFacVel = rampFactorVelSine(time, Tperiod);

		Real f, df;
		const Real bct     = t0 + t1 + t2 + t3;
		assert(bct>0);
		const Real shift   = std::floor((time-tStart)/bct);
		const Real tcoast  = tStart  + shift*bct;
		const Real tfreeze = tcoast  + t0;
		const Real tburst  = tfreeze + t1;
		const Real tswim   = tburst  + t2;
		//const Real phase   = (time<tfreeze) ?  shift   *0.5 + phaseShift
		//		 		  : (shift+1)*0.5 + phaseShift;
		const Real phase = 0.0;

		if (time<tcoast) {
			printf("NCCUCDC.\n");
			abort();
		} else if (time<tfreeze) {
			const Real d = (time-tcoast)/(tfreeze-tcoast);
			const std::pair<double, double> retVal = cubicHermite(1.0, lowestAmp, d);
			f = retVal.first;
			df = retVal.second;
			//f = 1 - 3*d*d + 2*d*d*d;
		} else if (time<tburst) {
			//f = 0.0;
			f = lowestAmp;
			df = 0.0;
		} else if (time<tswim) {
			const Real d = (time-tburst)/(tswim-tburst);
			const std::pair<double, double> retVal = cubicHermite(lowestAmp, 1.0, d);
			f = retVal.first;
			df = retVal.second;
			//f = 3*d*d - 2*d*d*d;
			//df = 6*(d - d*d)/(tswim-tburst);
		} else {
			f = 1.0;
			df = 0.0;
		}

		rX[0] = 0.0;
		rY[0] = rampFac*midlineBeC(rS[0], time, length, Tperiod, phase, f);
		for(int i=1;i<Nm;++i) {
			rY[i]=rampFac*midlineBeC(rS[i], time, length, Tperiod, phase, f);
			const Real dy = rY[i]-rY[i-1];
			const Real ds = rS[i] - rS[i-1];
			const Real dx = std::sqrt(ds*ds-dy*dy);
			rX[i] = rX[i-1] + dx;
		}


		vX[0] = 0.0; //rX[0] is constant
		vY[0] = rampFac*midlineVelBeC(rS[0],time,length,Tperiod, phase, f, df) +
			rampFacVel*midlineBeC(rS[0],time,length,Tperiod, phase, f);

		for(int i=1;i<Nm;++i) {
			vY[i]=rampFac*midlineVelBeC(rS[i],time,length,Tperiod, phase, f, df) +
				rampFacVel*midlineBeC(rS[i],time,length,Tperiod, phase, f);
			const Real dy  = rY[i]-rY[i-1];
			const Real dx  = rX[i]-rX[i-1];
			const Real dVy = vY[i]-vY[i-1];
			assert(dx>0); // has to be, otherwise y(s) is multiple valued for a given s
			vX[i] = vX[i-1] - dy/dx * dVy; // use ds^2 = dx^2 + dy^2 --> ddx = -dy/dx*ddy
		}

	}

	void _computeMidlineCoordinates(const Real time)
	{
		const Real rampFac = rampFactorSine(time, Tperiod);
		rX[0] = 0.0;
		rY[0] = rampFac*midline(rS[0], time, length, Tperiod, phaseShift);

		int hinge1Index = -1, hinge2Index = -1;

		for(int i=1;i<Nm;++i) {
			rY[i]=rampFac*midline(rS[i], time, length, Tperiod, phaseShift);
			const Real dy = rY[i]-rY[i-1];
			const Real ds = rS[i] - rS[i-1];
			Real dx = std::sqrt(ds*ds-dy*dy);

			// dx can be undef for s<0 and s>L points when wavelength>1. I dunno why we have these goddamn s<0 and s>L points
			if(not(dx>0) and not(waveLength==1.0)){
				dx = 0.001*length;
			}

			rX[i] = rX[i-1] + dx;

			if(rS[i]>=sHinge and hinge1Index<0) hinge1Index = i;
			if(rS[i]>=sHinge2 and hinge2Index<0) hinge2Index = i;
		}

		// Now do the second hinge section
		if(bDoubleHinge){

			//linearly decrease spring stiffness over 1 Tperiod, otherwise might get discontinuous theta2 at startup
			/*const bool kSpringTransition = std::floor(time-Tperiod) < 0;
			  const double kCurrent = not(kSpringTransition) ? kSpring : (kMaxSpring + time*(kSpring-kMaxSpring)/Tperiod);*/

			const double dt = time - this->timeNminus;
			this->timeNminus = time;
			this->runningTorque += this->torqueZsecMarkers * dt;

			if(time> (prevTransition+1)*0.01*Tperiod ){
				this->tOld = time;
				this->thetaOld = this->dTheta2;
				this->avgTorque = this->runningTorque/(0.01*Tperiod);
				this->runningTorque = 0.0;
				prevTransition++;
			}

			//Rigid until time period
			//if(time<Tperiod)
			if(1){
				this->dTheta2 = 0.0;
			}else{
				const double kCurrent = kSpring;
				const double cSpring = kCurrent;
				//const double thetaOld = this->dTheta2;
				const double thetaOld = this->thetaOld;
				const double tOld = this->tOld;
				const double torque = this->avgTorque;
				const double a1 = (thetaOld - torque/kCurrent)*exp(tOld*kCurrent/cSpring);
				//this->tOld = time;
				this->dTheta2 = torque/kCurrent + a1*exp(-time*kCurrent/cSpring);
				printf("time = %f, dTheta2 = %f, kSpring=%f, torque=%f\n", time, this->dTheta2, kCurrent, torque);
			}

			const double hinge1Loc[2] = {rX[hinge1Index], rY[hinge1Index]};
			const double hinge2Loc[2] = {rX[hinge2Index], rY[hinge2Index]};

			// angle of arm1 wrt main fish - imposed
			// Don't use analytical thetaExpression, since rampFactor not accounted-for in there
			const double currentTheta = std::atan( (hinge2Loc[1] - hinge1Loc[1]) / (hinge2Loc[0] - hinge1Loc[0]));

			for(int i=hinge2Index; i<Nm; ++i){
				const double dx = rX[i] - hinge2Loc[0];
				const double dy = rY[i] - hinge2Loc[1];
				const double localLength = std::sqrt(dx*dx + dy*dy);

				// angle of arm2 wrt main fish - from spring
				const double thetaHinge2 = currentTheta + dTheta2;

				rX[i] = hinge2Loc[0] + localLength*std::cos(thetaHinge2);
				rY[i] = hinge2Loc[1] + localLength*std::sin(thetaHinge2);
			}
		}
	}

	void _computeMidlineVelocities(const Real time)
	{
		const Real rampFac =    rampFactorSine(time, Tperiod);
		const Real rampFacVel = rampFactorVelSine(time, Tperiod);

		vX[0] = 0.0; //rX[0] is constant
		vY[0] = rampFac*midlineVel(rS[0],time,length,Tperiod, phaseShift) +
				rampFacVel*midline(rS[0],time,length,Tperiod, phaseShift);

		int indHinge2 = -1;
		for(int i=1;i<Nm;++i) {
			vY[i]=rampFac*midlineVel(rS[i],time,length,Tperiod, phaseShift) +
					rampFacVel*midline(rS[i],time,length,Tperiod, phaseShift);
			const Real dy  = rY[i]-rY[i-1];
			const Real dx  = rX[i]-rX[i-1];
			const Real dVy = vY[i]-vY[i-1];

			vX[i] = vX[i-1] - dy/dx * dVy; // use ds^2 = dx^2 + dy^2 --> ddx = -dy/dx*ddy
			if(waveLength==1.0){
				assert(dx>0); // has to be, otherwise y(s) is multiple valued for a given s
			}else{ // dx can be undef for s<0 and s>L points when wavelength>1. I dunno why we have these goddamn s<0 and s>L points
				if(not(dx>0))	vX[i] = 0.0;
			}

			if(indHinge2<0 and rS[i]>=sHinge2) indHinge2 = i;
		}

		if(bDoubleHinge){

			double dtHinge = time-oldTime;

			if(firstStep){
				for(int i=indHinge2; i<Nm; ++i){
					rXold[i] = rX[i];
					rYold[i] = rY[i];
				}
				firstStep = false;
				dtHinge = 1.0; // To avoid divide by zero at first step
			}

			for(int i=indHinge2; i<Nm; ++i){
				vX[i] = (rX[i] - rXold[i]) / dtHinge;
				vY[i] = (rY[i] - rYold[i]) / dtHinge;

				rXold[i] = rX[i];
				rYold[i] = rY[i];
			}
		}
		oldTime = time;

		/*FILE *temp;
		temp = fopen("vels.txt","a");
		for (int i=0; i<Nm; ++i){
			fprintf(temp,"%f\t", rS[i]);
		}

		fprintf(temp,"\n");
		for (int i=0; i<Nm; ++i){
			fprintf(temp,"%f\t", rX[i]);
		}

		fprintf(temp,"\n");
		for (int i=0; i<Nm; ++i){
			fprintf(temp,"%f\t", rY[i]);
		}

		fprintf(temp,"\n");
		for (int i=0; i<Nm; ++i){
			fprintf(temp,"%f\t", vX[i]);
		}

		fprintf(temp,"\n");
		for (int i=0; i<Nm; ++i){
			fprintf(temp,"%f\t", vY[i]);
		}

		fprintf(temp,"\n");
		fclose(temp);*/

	}

 public:
	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod,
		const Real phaseShift, const Real dx_ext, const Real _fac = 0.1212121212121212)
	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(false), tStart(1.0e12), bHinge(false), sHinge(0.0), AhingeTheta(0.0), ThingeTheta(0.0), hingePhi(0.0), quadraticAmplitude(false)
	{
	}

	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod,
		const Real phaseShift, const Real dx_ext, const string fburstpar, const Real _tStart, const Real _fac = 0.1212121212121212)
	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(true), tStart(_tStart), bHinge(false), sHinge(0.0), AhingeTheta(0.0), ThingeTheta(0.0), hingePhi(0.0), quadraticAmplitude(false)
	{
		//ifstream reader("burst_coast_carling_params.txt");
		ifstream reader(fburstpar.c_str());
		if (reader.is_open()) {
			reader >> t0;
			reader >> t1;
			reader >> t2;
			reader >> t3;
			reader >> lowestAmp;
			if(reader.eof()){
				cout << "Insufficient number of parameters provided for burstCoast" << endl;
				abort();
			}
			reader.close();
		} else {
			cout << "Could not open the correct 'params'.txt file" << endl;
			abort();
		}
	}

	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod,
			const Real phaseShift, const Real dx_ext, const Real _sHinge,
			const double _Ahinge, const double _phiHinge, const Real _fac = 0.1212121212121212):
		FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(false),tStart(1e12), bHinge(true), sHinge(_sHinge), AhingeTheta(M_PI*_Ahinge/180.0), ThingeTheta(Tperiod), hingePhi(_phiHinge/360.0), quadraticAmplitude(true)
    {
    }

	// Had to get rid of default value of _fac=0.1212121212121212. Was creating confusion for compiler with plain fish constructor 
	CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod,
			const Real phaseShift, const Real dx_ext, const Real _sHinge,
			const Real _fac):
		FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(false),tStart(1e12), bHinge(true), sHinge(_sHinge), ThingeTheta(Tperiod), quadraticAmplitude(true)
	{
		// Now, read the optimization params (Ahinge, _phiHinge, tail size) from the params file
		{
			Real _Ahinge, _phiHinge;
			//ifstream reader("burst_coast_carling_params.txt");
			ifstream reader("hingedParams.txt");
			if (reader.is_open()) {
				reader >> _Ahinge;
				reader >> _phiHinge;
				reader >> finSize;
				reader >> waveLength;
				reader >> sHinge2;
				reader >> kSpring;
				if(reader.eof()){
					cout << "Insufficient number of parameters provided for hingedFin" << endl; fflush(NULL); abort();
				}
				reader.close();
			} else {
				cout << "Could not open the correct 'params'.txt file" << endl; fflush(NULL);
				abort();
			}

			AhingeTheta = M_PI*_Ahinge/180.0;
			hingePhi = _phiHinge/360.0;
			sHinge2 *= length;
			// UnRescaling: to avoid CMA trouble
			kSpring *= 1.0e-4;
		}

		// FinSize has now been updated with value read from text file. Recompute heights to over-write with updated values
		_computeWidthsHeights();

	}

	void computeMidline(const Real time) override
	{

		if (bBurst && time>=tStart) {
			_computeMidlineCoordinatesBeC(time);
		} else {
			_computeMidlineCoordinates(time);
			_computeMidlineVelocities(time);
		}
		_computeMidlineNormals();
		#if 0
		#warning USED MPI COMM WORLD
				// we dump the profile
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD,&rank);
				if (rank!=0) return;
				FILE * f = fopen("fish_profile","w");
				for(int i=0;i<Nm;++i)
					fprintf(f,"%d %g %g %g %g %g %g %g %g %g %g %g\n",
							i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i],
							vNorX[i],vNorY[i],width[i],height[i]);
				fclose(f);
				printf("Dumped midline\n");
		#endif
	}
};

class CurvatureDefinedFishData : public FishMidlineData
{
protected:
	Real * const rK;
	Real * const vK;
	Real * const rC;
	Real * const vC;
	Real * const rB;
	Real * const vB;
	Real * const rA;
	Real * const vA;
	Real controlFac, valPID;
	Real controlVel, velPID;
public:

	CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext)
	: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext),
	  rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)),
		rA(_alloc(Nm)),vA(_alloc(Nm)), rB(_alloc(Nm)),vB(_alloc(Nm)), controlFac(-1), controlVel(0), velPID(0),valPID(0)
	{
	}

	void _correctTrajectory(const Real dtheta, const Real vtheta, const Real time, Real dt) override
	{
		velPID = vtheta;
		valPID = dtheta;
		
		dt = std::max(2.2e-16,dt);
		std::array<Real,6> tmp_curv = std::array<Real,6>();
		for (int i=0; i<tmp_curv.size(); ++i) {tmp_curv[i] = dtheta;}
		//adjustScheduler.transition(time,time,time+2*dt,tmp_curv, true);
		adjustScheduler.transition(time, time-2*dt, time+2*dt, tmp_curv, true);
		
	}

	void _correctAmplitude(Real dAmp, Real vAmp, const Real time, const Real dt) override
	{
		assert(dAmp>0 && dAmp<2); //buhu
		if(dAmp<=0) {
			dAmp=0; 
			vAmp=0;
		}
		controlFac = dAmp;
		controlVel = vAmp;
		//const Real rampUp = time<Tperiod ? time/Tperiod : 1; //TODO actually should be cubic spline!
		//const Real fac = dAmp*rampUp/length; //curvature is 1/length
		//const std::array<Real ,6> curvature_values = {
		//	fac*0.82014, fac*1.46515, fac*2.57136, fac*3.75425, fac*5.09147, fac*5.70449
		//};
		//curvScheduler.transition(time,time,time+2*dt, curvature_values, true);
		//curvScheduler.transition(time, time-dt, time+dt, curvature_values);
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
			printf("Turning by %g at time %g with period %g.\n", input[0], time, l_tnext);
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

	void computeMidline(const Real time) override
	{
		const Real _1oL = 1./length;
		const Real _1oT = 1./l_Tp;
		const std::array<Real ,6> curvature_points = {
				0, .15*length, .4*length, .65*length, .9*length, length
		};
		const std::array<Real,7> baseline_points = {-.5,-.25,0,.25,.5,.75,1.};
		const std::array<Real ,6> curvature_values = {
				0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
				3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
		};
		const std::array<Real,6> curvature_zeros = std::array<Real, 6>();
		curvScheduler.transition(time,0,Tperiod,curvature_zeros,curvature_values);

		// query the schedulers for current values
		curvScheduler.gimmeValues(  time, 							curvature_points, Nm, rS, rC, vC);
		baseScheduler.gimmeValues(  time, l_Tp, length, baseline_points, 	Nm, rS, rB, vB);
		adjustScheduler.gimmeValues(time, 							curvature_points, Nm, rS, rA, vA);
		if(controlFac>0) {
			const Real _vA = velPID, _rA = valPID;
			// construct the curvature
			for(unsigned int i=0; i<Nm; i++) {
				const Real darg = 2.*M_PI* _1oT;
				const Real arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL/waveLength) + M_PI*phaseShift;
				rK[i] =   rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlFac;
				vK[i] =   vC[i]*(std::sin(arg)     +rB[i]+_rA)*controlFac
					+ rC[i]*(std::cos(arg)*darg+vB[i]+_vA)*controlFac
					+ rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlVel;
			}
		} else {
			// construct the curvature
			for(unsigned int i=0; i<Nm; i++) {
				const Real darg = 2.*M_PI* _1oT;
				const Real arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL/waveLength) + M_PI*phaseShift;
				rK[i] =   rC[i]*(std::sin(arg)      + rB[i] + rA[i]);
				vK[i] =   vC[i]*(std::sin(arg)      + rB[i] + rA[i])
					+ rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]);
			}
		}


		//printf("%g %g %g %g\n", rB[12], rB[425], rB[838], rB[1238]);



		#if 0
			{ // we dump the profile points
				FILE * f = fopen("stefan.dat","a");
				std::array<Real, 6> curv,base;
				curvScheduler.ParameterScheduler<6>::gimmeValues(time, curv);
				baseScheduler.ParameterScheduler<6>::gimmeValues(time, base);
				fprintf(f,"%9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e\n",
						time,curv[0],curv[1],curv[2],curv[3],curv[4],curv[5],base[0],base[1],base[2],base[3],base[4],base[5]);
				fclose(f);
			}
		#endif
		// solve frenet to compute midline parameters
		IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
		#if 0
			#warning USED MPI COMM WORLD
			{
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD,&rank);
				if (rank!=0) return;
				FILE * f = fopen("stefan_profile","w");
				for(int i=0;i<Nm;++i)
					fprintf(f,"%d %g %g %g %g %g %g %g %g %g\n",
						i,rS[i],rX[i],rY[i],vX[i],vY[i],
						vNorX[i],vNorY[i],width[i],height[i]);
				fclose(f);
			}
		#endif

	}
};

struct VolumeSegment_OBB
{
	std::pair<int, int> s_range;
	Real normalI[3]; // should be normalized and >=0
	Real normalJ[3];
	Real normalK[3];
	Real w[3]; // halfwidth
	Real c[3]; // center

	#if 1==0
	VolumeSegment_OBB(std::pair<int, int> s_range, const Real bbox[3][2])
	: s_range(s_range)
	{
		normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
		normalI[0]=normalJ[1]=normalK[2]=1.0;
		for(int i=0;i<3;++i) {
			w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
			c[i] = bbox[i][0] + w[i];
			assert(w[i]>0);
		}

	}
	#endif

	VolumeSegment_OBB() { }

	void prepare(std::pair<int, int> _s_range, const Real bbox[3][2])
	{
		normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
		normalI[0]=normalJ[1]=normalK[2]=1.0;
		s_range.first = _s_range.first;
		s_range.second = _s_range.second;
		for(int i=0;i<3;++i) {
			w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
			c[i] = bbox[i][0] + w[i];
			assert(w[i]>0);
		}
	}

	void normalizeNormals()
	{
		const Real magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
		const Real magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
		const Real magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
		assert(magI > std::numeric_limits<Real>::epsilon());
		assert(magJ > std::numeric_limits<Real>::epsilon());
		assert(magK > std::numeric_limits<Real>::epsilon());
		const Real invMagI = 1.0/magI;
		const Real invMagJ = 1.0/magJ;
		const Real invMagK = 1.0/magK;

		for(int i=0;i<3;++i) {
			// also take absolute value since thats what we need when doing intersection checks later
			normalI[i]=std::fabs(normalI[i])*invMagI;
			normalJ[i]=std::fabs(normalJ[i])*invMagJ;
			normalK[i]=std::fabs(normalK[i])*invMagK;
		}
	}

	void changeToComputationalFrame(const Real position[3], const Real quaternion[4])
	{
		// we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
		const Real w = quaternion[0];
		const Real x = quaternion[1];
		const Real y = quaternion[2];
		const Real z = quaternion[3];
		const Real Rmatrix3D[3][3] = {
				{1.-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
				{2*(x*y+z*w),    1.-2*(x*x+z*z),  2*(y*z-x*w)},
				{2*(x*z-y*w),    2*(y*z+x*w),    1.-2*(x*x+y*y)}
		};
		const Real p[3] = {c[0],c[1],c[2]};
		const Real nx[3] = {normalI[0],normalI[1],normalI[2]};
		const Real ny[3] = {normalJ[0],normalJ[1],normalJ[2]};
		const Real nz[3] = {normalK[0],normalK[1],normalK[2]};
		for(int i=0;i<3;++i) {
			c[i] = Rmatrix3D[i][0]*p[0] + Rmatrix3D[i][1]*p[1] + Rmatrix3D[i][2]*p[2];
			normalI[i] = Rmatrix3D[i][0]*nx[0] + Rmatrix3D[i][1]*nx[1] + Rmatrix3D[i][2]*nx[2];
			normalJ[i] = Rmatrix3D[i][0]*ny[0] + Rmatrix3D[i][1]*ny[1] + Rmatrix3D[i][2]*ny[2];
			normalK[i] = Rmatrix3D[i][0]*nz[0] + Rmatrix3D[i][1]*nz[1] + Rmatrix3D[i][2]*nz[2];
		}
		c[0] +=position[0];
		c[1] +=position[1];
		c[2] +=position[2];
		normalizeNormals();
	}

	bool isIntersectingWithAABB(const Real start[3],const Real end[3], const Real h_gridpt) const
	{
		//start and end are two diagonally opposed corners of grid block
		/*const Real AABB_w[3] = { //half block width + safe distance
		  0.5*(end[0] - start[0]) + 2.0*safe_distance,
		  0.5*(end[1] - start[1]) + 2.0*safe_distance,
		  0.5*(end[2] - start[2]) + 2.0*safe_distance
		  };*/
		const Real AABB_w[3] = { //half block width + 1 point on either side (Keep extra point, since cell is empty between center and face)
			0.5*(end[0] - start[0] + h_gridpt) + h_gridpt, // Need (end-start+h) since 'start' and 'end' correspond to cell centers, not the faces
			0.5*(end[1] - start[1] + h_gridpt) + h_gridpt,
			0.5*(end[2] - start[2] + h_gridpt) + h_gridpt
		};
		/*const Real AABB_c[3] = { //block center
		  start[0] + AABB_w[0] - safe_distance,
		  start[1] + AABB_w[1] - safe_distance,
		  start[2] + AABB_w[2] - safe_distance
		  };*/
		const Real AABB_c[3] = { //block center
			0.5*(start[0] + end[0]),
			0.5*(start[1] + end[1]),
			0.5*(start[2] + end[2])
		};

		assert(AABB_w[0]>0);
		assert(AABB_w[1]>0);
		assert(AABB_w[2]>0);
		bool intersects = true;

		const Real r1 = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
		if (not ((c[0]-r1 <= AABB_c[0] + AABB_w[0]) && (c[0]+r1 >= AABB_c[0] - AABB_w[0])))
			return false;

		const Real r2 = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
		if (not ((c[1]-r2 <= AABB_c[1] + AABB_w[1]) && (c[1]+r2 >= AABB_c[1] - AABB_w[1])))
			return false;

		const Real r3 = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
		if (not ((c[2]-r3 <= AABB_c[2] + AABB_w[2]) && (c[2]+r3 >= AABB_c[2] - AABB_w[2])))
			return false;

		const Real r4 = AABB_w[0]*normalI[0] + AABB_w[1]*normalI[1] + AABB_w[2]*normalI[2];
		//const Real r4 = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
		if (not ((AABB_c[0]-r4 <= c[0] + w[0]) && (AABB_c[0]+r4 >= c[0] - w[0])))
			return false;

		const Real r5 = AABB_w[0]*normalJ[0] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalJ[2];
		//const Real r5 = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
		if (not ((AABB_c[1]-r5 <= c[1] + w[1]) && (AABB_c[1]+r5 >= c[1] - w[1])))
			return false;

		const Real r6 = AABB_w[0]*normalK[0] + AABB_w[1]*normalK[1] + AABB_w[2]*normalK[2];
		//const Real r6 = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
		if (not ((AABB_c[2]-r6 <= c[2] + w[2]) && (AABB_c[2]+r6 >= c[2] - w[2])))
			return false;

		return true;
	}
};

struct PutFishOnBlocks
{
	const FishMidlineData * cfish;
	Real position[3];
	Real quaternion[4];
	Real Rmatrix3D[3][3];

	PutFishOnBlocks(const FishMidlineData* const cfish, const Real pos[3], const Real quat[4]):
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

	void computeRotationMatrix()
	{
		const Real w = quaternion[0];
		const Real x = quaternion[1];
		const Real y = quaternion[2];
		const Real z = quaternion[3];
		const Real R[3][3] = {
				{1-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
				{2*(x*y+z*w),    1-2*(x*x+z*z),  2*(y*z-x*w)},
				{2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x*x+y*y)}
		};
		memcpy(Rmatrix3D, R, sizeof(R));
	}

	inline int find_closest_dist(const int s, const int dir, const Real x[3], Real & oldDistSq) const
	{
		if((s+dir)<cfish->iFishStart or (s+dir)>cfish->iFishEnd)
			return s;

		const Real newDistSq = (x[0]-cfish->rX[s+dir])*(x[0]-cfish->rX[s+dir])
											   + (x[1]-cfish->rY[s+dir])*(x[1]-cfish->rY[s+dir])
																											   + (x[2])*(x[2]);

		if(oldDistSq<=newDistSq) {
			return s;
		} else {
			oldDistSq = newDistSq;
			return s+dir;
		}
	}

	void changeVelocityToComputationalFrame(Real x[3]) const
	{
		const Real p[3] = {x[0],x[1],x[2]};

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

	Real getSmallerDistToMidline(const int start_s, const Real x[3], int & final_s) const
	{
		Real relX[3] = {x[0],x[1],x[2]};
		changeFromComputationalFrame(relX);

		const Real curDistSq =  std::pow(relX[0]-cfish->rX[start_s],2)
		+ std::pow(relX[1]-cfish->rY[start_s],2)
		+ std::pow(relX[2],2);

		Real distSq;

		distSq = curDistSq; // check right
		const int sRight = find_closest_dist(start_s, +1, relX, distSq);

		distSq = curDistSq; // check left
		const int sLeft = find_closest_dist(start_s, -1, relX, distSq);

		if(sRight==start_s and sLeft==start_s) {
			final_s = start_s;
			return distSq;
		}

		/*if (not (sRight==start_s or sLeft==start_s)){
			FILE * pFile;
			pFile = fopen("problematique.txt","w");
			fprintf(pFile, "curDistSq = %18.16f, distSq = %18.16f\n", curDistSq, distSq);
			fprintf(pFile, "targetX= %18.16f, targetY= %18.16f\n", relX[0], relX[1]);
			fprintf(pFile, "sLeft=%d, start_s=%d, sRight=%d\n", sLeft, start_s, sRight);
			for (int i=0; i<cfish.Nm; i++){
				fprintf(pFile, "%18.16f \t %18.16f\n", cfish.rX[i], cfish.rY[i]);
			}
			fclose(pFile);
		}*/

		// With a hinged tail, the following overzealous assert will catch, outputted problematique and confirmed in Matlab
		//assert(sRight==start_s or sLeft==start_s);

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
			Real * const chi = &(defblock->chi[0][0][0]);
			Real * const udef = &(defblock->udef[0][0][0][0]);

			static const int n = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
			for(int i=0; i<n; i++) {
				chi[i]=-1;
				//udef[3*i+0]=0;
				//udef[3*i+1]=0;
				//udef[3*i+2]=0;
			}
		}

		// construct the shape (P2M with min(distance) as kernel) onto defblocks
		for(int i=0;i<vSegments.size();++i) {
			//iterate over segments contained in the vSegm intersecting this block:
			const int firstSegm = std::max(vSegments[i].s_range.first,cfish->iFishStart);
			const int lastSegm = std::min(vSegments[i].s_range.second, cfish->iFishEnd);
			for(int ss=firstSegm; ss<=lastSegm; ++ss) {
				assert(ss>=cfish->iFishStart && ss<=cfish->iFishEnd);

				// fill chi
				// assume width is major axis, else correction:
				const Real offset = cfish->height[ss] > cfish->width[ss] ? 0.5*M_PI : 0.0;
				const Real ell_a = (Real)std::max(cfish->height[ss],cfish->width[ss]);
				//                    const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : 0.25*info.h[0]/ell_a;
				// maximum distance between two points is 2*ell_a * sin(theta). set this distance to dx/2 -->
				const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : std::abs(std::asin(0.5*(0.25*info.h_gridpoint)/ell_a));

				const int Ntheta = (int)std::ceil(2.0*M_PI/dtheta_target) + 1;
				const Real dtheta = 2.0*M_PI/((Real)Ntheta);

				for(int tt=0;tt<Ntheta;++tt) {
					const Real theta = tt*dtheta + offset;
					// create a surface point
					// special treatment of tail (width = 0 --> no ellipse, just line)
					const Real hght = cfish->width[ss] == 0 ? cfish->height[ss]*(2*tt/((Real)Ntheta-1) - 1)
							: cfish->height[ss]*std::sin(theta);
					Real myP[3] = {
							(cfish->rX[ss] + cfish->width[ss]*std::cos(theta)*cfish->norX[ss]),
							(cfish->rY[ss] + cfish->width[ss]*std::cos(theta)*cfish->norY[ss]),
							hght
					};
					changeToComputationalFrame(myP);

					const int iap[3] = {
							(int)std::floor((myP[0]-org[0])*invh),
							(int)std::floor((myP[1]-org[1])*invh),
							(int)std::floor((myP[2]-org[2])*invh)
					};


					// Don't try to populate outside bounds
                                        const int startMarker[3] = {
                                                        std::max(0,iap[0]),
                                                        std::max(0,iap[1]),
                                                        std::max(0,iap[2])
                                        };
                                        const int endMarker[3] = {
                                                        std::min(iap[0]+2, FluidBlock::sizeX - 0), //+2 instead of +1, since will do (< than)
                                                        std::min(iap[1]+2, FluidBlock::sizeY - 0), //weird: if I don't put -0, FluidBlock not recognized
                                                        std::min(iap[2]+2, FluidBlock::sizeZ - 0)
                                        };

					//Prep section, finalize later towards end (all ellipse encompassing points grabbed if beyond hinge2)
					if(ss>=defblock->hinge2Index){
						// Populate in the 'box' containing the current point. Low corners are iap, high corners are iap+1
						for(int sz=startMarker[2]; sz<endMarker[2]; ++sz){
							for(int sy=startMarker[1]; sy<endMarker[1]; ++sy){
								for(int sx=startMarker[0]; sx<endMarker[0]; ++sx){
									defblock->sectionMarker[sz][sy][sx] = 1;
								}
							}
						}
					}


					// support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
					const int start[3] = {
							std::max(-1, 0 - iap[0] ),
							std::max(-1, 0 - iap[1] ),
							std::max(-1, 0 - iap[2] )
					};
					const int end[3] = {
							std::min(+3, FluidBlock::sizeX - iap[0]),
							std::min(+3, FluidBlock::sizeY - iap[1]),
							std::min(+3, FluidBlock::sizeZ - iap[2])
					};
					const Real myP_distToMidlineSq = cfish->width[ss] == 0 ? std::pow(hght,2) :
							(Real)(std::pow(cfish->width[ss]*std::cos(theta),2) + std::pow(cfish->height[ss]*std::sin(theta),2));

					if(myP_distToMidlineSq<std::numeric_limits<Real>::epsilon()) {
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
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);

							changeFromComputationalFrame(p);
							const Real distPlanar = std::sqrt( std::pow(p[0]-cfish->rX[closest_s],2) +
									std::pow(p[1]-cfish->rY[closest_s],2) );
							const Real distHeight = std::abs(p[2]);
							const Real sign = (distPlanar > cfish->width[closest_s] or distHeight > cfish->height[closest_s]) ? -1.0 : 1.0;

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
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);
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
				assert(ss>=0 && ss<=cfish->Nm-1);
				// P2M udef of a slice at this s
				const Real myWidth =  (ss < cfish->iFishStart ? cfish->width[ cfish->iFishStart]
												  	: (ss > cfish->iFishEnd   ? cfish->width[ cfish->iFishEnd]
																											: cfish->width[ss]));
				const Real myHeight = (ss < cfish->iFishStart ? cfish->height[cfish->iFishStart]
													  : (ss > cfish->iFishEnd   ? cfish->height[cfish->iFishEnd]
																											: cfish->height[ss]));
				const Real ds_defGrid = info.h_gridpoint;
				// towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice (and so we can floor!)
				const Real extension = NPPEXT*info.h_gridpoint; //G tmp changed back to 2

				const int Nh = std::floor( (myHeight+extension)/ds_defGrid );

				for(int ih=-Nh;ih<=Nh; ++ih) {
					const Real offsetH = ih*ds_defGrid;
					// add an extra extension when width == 0 (to deal with large curvatures near head and/or tail):
					const Real currentWidth = myWidth== 0 ? extension : myWidth * std::sqrt(1 - std::pow(offsetH/(myHeight+extension),2));
					const Real actualWidth = (cfish->height[ss] == 0 or std::abs(offsetH)>=cfish->height[ss]) ? 0.0
							: cfish->width[ss] * std::sqrt(1 - std::pow(offsetH/cfish->height[ss],2));
					const int Nw = std::floor( (currentWidth+extension)/ds_defGrid); // add extension here to make sure we have it in each direction

					for(int iw=-Nw;iw<=Nw; ++iw) {
						const Real offsetW = iw*ds_defGrid;
						Real xp[3] = {
								(cfish->rX[ss] + offsetW*cfish->norX[ss]),
								(cfish->rY[ss] + offsetW*cfish->norY[ss]),
								offsetH
						};
						changeToComputationalFrame(xp);
						xp[0] = (xp[0]-org[0])*invh;
						xp[1] = (xp[1]-org[1])*invh;
						xp[2] = (xp[2]-org[2])*invh;
						Real udef[3] = {
								(cfish->vX[ss] + offsetW*cfish->vNorX[ss]),
								(cfish->vY[ss] + offsetW*cfish->vNorY[ss]),
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
								std::min(+2, FluidBlock::sizeX - iap[0]),
								std::min(+2, FluidBlock::sizeY - iap[1]),
								std::min(+2, FluidBlock::sizeZ - iap[2])
						};

						Real wghts[3][2];
						for(int c=0;c<3;++c) {
							const Real t[2] = {
									std::fabs((Real)xp[c] - (ap[c]+0)),
									std::fabs((Real)xp[c] - (ap[c]+1))
							};
							wghts[c][0] = 1.0 - t[0];
							wghts[c][1] = 1.0 - t[1];
						}

						const bool isInside =  (std::fabs(offsetW) < actualWidth)
																&& (std::fabs(offsetH) < cfish->height[ss]);
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
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);
							defblock->udef[idx[2]][idx[1]][idx[0]][0] += wxwywz*udef[0];
							defblock->udef[idx[2]][idx[1]][idx[0]][1] += wxwywz*udef[1];
							defblock->udef[idx[2]][idx[1]][idx[0]][2] += wxwywz*udef[2];
							b(idx[0],idx[1],idx[2]).tmpU += wxwywz;
							// set sign for all interior points
							if( (std::fabs(defblock->chi[idx[2]][idx[1]][idx[0]] + 1) <
											5*std::numeric_limits<Real>::epsilon()) && isInside)
								defblock->chi[idx[2]][idx[1]][idx[0]] = 1.0;
						}
							}
						}
					}
				}
			}
		}

		// finalize signed distance function in tmpU
		{
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++) {
				const Real normfac = b(ix,iy,iz).tmpU > numeric_limits<Real>::epsilon() ? b(ix,iy,iz).tmpU : 1.;
				defblock->udef[iz][iy][ix][0] /= normfac;
				defblock->udef[iz][iy][ix][1] /= normfac;
				defblock->udef[iz][iy][ix][2] /= normfac;
				// change from signed squared distance function to normal sdf
				b(ix,iy,iz).tmpU = defblock->chi[iz][iy][ix] > 0 ?
						sqrt( defblock->chi[iz][iy][ix]) :
						-sqrt(-defblock->chi[iz][iy][ix]);
				b(ix,iy,iz).tmpV = defblock->udef[iz][iy][ix][0];
				b(ix,iy,iz).tmpW = defblock->udef[iz][iy][ix][1];

				// All points that are not chi=0 in the targeted section, are captured here. When we loop through SurfaceBlocks for computing torque, the extraneous points captured here will be left out, so hakunamatata.
				defblock->sectionMarker[iz][iy][ix] *= std::abs(defblock->chi[iz][iy][ix]) > 0 ? 1 : 0;
			}
		}
	}
};

struct PutFishOnBlocks_Finalize : public GenericLabOperator
{
	Real t;
	int stencil_start[3], stencil_end[3];
	array<Real,4>* const momenta;
	surfaceBlocks* const surface;
	std::map<int,ObstacleBlock*>* const obstacleBlocks;
	//PutFishOnBlocks_Finalize finalize(obstacleBlocks,dataPerThread[tid],tmp,blockID);

	PutFishOnBlocks_Finalize(map<int,ObstacleBlock*>* const obstacleBlocks, //to write chi
			surfaceBlocks* const surface, 					//to write gradChi
			array<Real,4>* const momenta)     			//compute CM
	: t(0), momenta(momenta),surface(surface), obstacleBlocks(obstacleBlocks)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 5);
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}
	#if 0
	inline Real sign(const Real& val) const {
		return (0. < val) - (val < 0.);
	}

	template <typename Lab, typename BlockType>
	void operator()(Lab& lab, const BlockInfo& info, BlockType& b)
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;

		const Real eps  = std::numeric_limits<Real>::epsilon();
		const Real h    = info.h_gridpoint;
		const Real fac1 = .5/h;
		const Real fac2 = .5/(std::sqrt(2.)*h);
		const Real fac[9] = {fac1,fac1,fac1,fac2,fac2,fac2,fac2,fac2,fac2};

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU > +2*h || lab(ix,iy,iz).tmpU < -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				(*momenta)[0] += H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				continue;
			}
			const Real Up[9] = {
					lab(ix+1,iy  ,iz  ).tmpU, lab(ix  ,iy+1,iz  ).tmpU, lab(ix  ,iy  ,iz+1).tmpU,
					lab(ix  ,iy+1,iz+1).tmpU, lab(ix+1,iy  ,iz+1).tmpU, lab(ix+1,iy+1,iz  ).tmpU,
					lab(ix  ,iy+1,iz-1).tmpU, lab(ix+1,iy  ,iz-1).tmpU, lab(ix+1,iy-1,iz  ).tmpU
			};
			const Real Um[9] = {
					lab(ix-1,iy  ,iz  ).tmpU, lab(ix  ,iy-1,iz  ).tmpU, lab(ix  ,iy  ,iz-1).tmpU,
					lab(ix  ,iy-1,iz-1).tmpU, lab(ix-1,iy  ,iz-1).tmpU, lab(ix-1,iy-1,iz  ).tmpU,
					lab(ix  ,iy-1,iz+1).tmpU, lab(ix-1,iy  ,iz+1).tmpU, lab(ix-1,iy+1,iz  ).tmpU
			};
			Real Ip[9],Im[9];
			for (int i=0; i<9; i++) {
				Ip[i] = std::max(Up[i],(Real)0.);
				Im[i] = std::max(Um[i],(Real)0.);
			}
			Real gradU[9], gradI[9];
			for (int i=0; i<9; i++) {
				gradU[i] = fac[i]*(Up[i]-Um[i]);
				gradI[i] = fac[i]*(Ip[i]-Im[i]);
			}
			Real Hp[3], Hm[3]; //only x y and z, the direction where i need gradChi
			for(int i=0; i<3; i++) {
				Hp[i] = (Up[i]> h) ? h : (
								(Up[i]<-h) ? 0 :
														.5*h + (Up[i] - fac1*sign(Up[i])*Up[i]*Up[i]) );
				Hm[i] = (Um[i]> h) ? h : (
								(Um[i]<-h) ? 0 :
														.5*h + (Um[i] - fac1*sign(Um[i])*Um[i]*Um[i]) );
			}

			const Real gradH[3] = {.5*(Hp[0]-Hm[0]), .5*(Hp[1]-Hm[1]), .5*(Hp[2]-Hm[2])};
			Real gradUU[3], gradUI[3], gradUH[3];
			for (int i=0; i<3; i++) {
				gradUU[i] = gradU[i]*gradU[i] + gradU[i+3]*gradU[i+3] + gradU[i+6]*gradU[i+6];
				gradUI[i] = gradU[i]*gradI[i] + gradU[i+3]*gradI[i+3] + gradU[i+6]*gradI[i+6];
				gradUH[i] = gradU[i]*gradH[i];
			}
			for (int i=0; i<3; i++)  gradUU[i] = max(gradUU[i], eps);
			const Real FDH = 1/3. *(gradUI[0]/gradUU[0]+gradUI[1]/gradUU[1]+gradUI[2]/gradUU[2]);

			const Real gradUSq = gradU[0]*gradI[0] + gradU[1]*gradI[1] + gradU[2]*gradI[2];
			const Real FDD = gradUSq<eps?0:(gradUH[0]+gradUH[1]+gradUH[2])*h/gradUSq; // == delta * h^3

			if (FDD>eps) {
				const Real dchidx = -FDD*gradU[0];
				const Real dchidy = -FDD*gradU[1];
				const Real dchidz = -FDD*gradU[2];
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, FDD);
			}
			#ifndef NDEBUG
			//		if(FDH<0 || FDH>1) printf("invalid H?: %9.9e %9.9e %9.9e: %9.9e\n",x,y,z,FDH);
			#endif
			(*momenta)[0] += FDH;
			(*momenta)[1] += p[0]*FDH;
			(*momenta)[2] += p[1]*FDH;
			(*momenta)[3] += p[2]*FDH;
			defblock->chi[iz][iy][ix] = FDH;
			b(ix,iy,iz).chi = std::max(FDH, b(ix,iy,iz).chi);
		}
	}
	#else
	template <typename Lab, typename BlockType>
	void operator()(Lab& lab, const BlockInfo& info, BlockType& b)
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;
		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real h = info.h_gridpoint;
		const Real inv2h = .5/h;
		const Real fac1 = 0.5*h*h;
		const Real fac2 = h*h*h;

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU > +2*h || lab(ix,iy,iz).tmpU < -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				(*momenta)[0] += H;
				continue;
			}

			const Real distPx = lab(ix+1,iy,iz).tmpU;
			const Real distMx = lab(ix-1,iy,iz).tmpU;
			const Real distPy = lab(ix,iy+1,iz).tmpU;
			const Real distMy = lab(ix,iy-1,iz).tmpU;
			const Real distPz = lab(ix,iy,iz+1).tmpU;
			const Real distMz = lab(ix,iy,iz-1).tmpU;
			// gradU
			const Real gradUX = inv2h*(distPx - distMx);
			const Real gradUY = inv2h*(distPy - distMy);
			const Real gradUZ = inv2h*(distPz - distMz);
			const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ + eps;

			/*
			if (gradUSq < eps) {
				b(ix,iy,iz).chi = std::max((Real)0, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = 0;
				continue;
			}
			*/
			const Real IplusX = distPx < 0 ? 0 : distPx;
			const Real IminuX = distMx < 0 ? 0 : distMx;
			const Real IplusY = distPy < 0 ? 0 : distPy;
			const Real IminuY = distMy < 0 ? 0 : distMy;
			const Real IplusZ = distPz < 0 ? 0 : distPz;
			const Real IminuZ = distMz < 0 ? 0 : distMz;
			const Real HplusX = distPx == 0 ? 0.5 : (distPx < 0 ? 0 : 1);
			const Real HminuX = distMx == 0 ? 0.5 : (distMx < 0 ? 0 : 1);
			const Real HplusY = distPy == 0 ? 0.5 : (distPy < 0 ? 0 : 1);
			const Real HminuY = distMy == 0 ? 0.5 : (distMy < 0 ? 0 : 1);
			const Real HplusZ = distPz == 0 ? 0.5 : (distPz < 0 ? 0 : 1);
			const Real HminuZ = distMz == 0 ? 0.5 : (distMz < 0 ? 0 : 1);

			// gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
			const Real gradIX = inv2h*(IplusX - IminuX);
			const Real gradIY = inv2h*(IplusY - IminuY);
			const Real gradIZ = inv2h*(IplusZ - IminuZ);
			const Real gradHX = (HplusX - HminuX);
			const Real gradHY = (HplusY - HminuY);
			const Real gradHZ = (HplusZ - HminuZ);
			const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;
			const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;
			const Real Delta = numD/gradUSq; //h^3 * Delta
			const Real H     = numH/gradUSq;

			if (Delta>1e-6) {
				const Real _Delta =  Delta*fac1;
				const Real dchidx = -Delta*gradUX*fac1;
				const Real dchidy = -Delta*gradUY*fac1;
				const Real dchidz = -Delta*gradUZ*fac1;
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, _Delta);
			}
			(*momenta)[0] += H;
			(*momenta)[1] += p[0]*H;
			(*momenta)[2] += p[1]*H;
			(*momenta)[3] += p[2]*H;
			defblock->chi[iz][iy][ix] = H;
			b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
		}
	}
	#endif
};

#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
