//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_FishLibrary__
#define __IncompressibleFluids3D__IF3D_FishLibrary__

#include <cmath>
#include <array>
#include "IF2D_Frenet.h"

const int NPPSEG = 50.; //was 100
const int NPPEXT = 3; //was 3
const int TGTPPB = 4.; //was 2 i think
const int TSTART = 2.;

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
			savestream << parameters_t0[i] << "\t" << parameters_t1[i] << "\t" << dparameters_t0[i] << std::endl;
		savestream.close();
	}

	void restart(std::string filename)
	{
		std::ifstream restartstream;
		restartstream.open(filename+".txt");

		restartstream >> t0 >> t1;
		for(int i=0;i<Npoints;++i)
			restartstream >> parameters_t0[i] >> parameters_t1[i] >> dparameters_t0[i];
		restartstream.close();
	}

	ParameterScheduler()
	{
		t0=-1;
		parameters_t0 = std::array<Real, Npoints>();
		dparameters_t0 = std::array<Real, Npoints>();
	}

	void transition(const Real t, const Real tstart, const Real tend,
			const std::array<Real, Npoints> parameters_tend,
			const bool UseCurrentDerivative = false)
	{
		if(t<tstart or t>tend) return; // this transition is out of scope
		if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

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
	void gimmeValues(const Real t, const Real Twave, const Real Length, const std::array<Real, Npoints> & positions,
			const int Nfine, const Real* const positions_fine, Real* const parameters_fine, Real* const dparameters_fine)
	{
		const Real _1oL = 1./Length;
		const Real _1oT = 1./Twave;
		// the fish goes through (as function of t and s) a wave function that describes the curvature
		for(int i=0;i<Nfine;++i) {
			const Real c = positions_fine[i]*_1oL - (t - this->t0)*_1oT; //traveling wave coord
			bool bCheck = true;

			if (c < positions[0]) { // Are you before latest wave node?
				IF2D_Interpolation1D::cubicInterpolation(c, positions[0], c,
						this->parameters_t0[0], this->parameters_t0[0],
						parameters_fine[i], dparameters_fine[i]);
				bCheck = false;
			}
			else if (c > positions[Npoints-1]) {// Are you after oldest wave node?
				IF2D_Interpolation1D::cubicInterpolation(positions[Npoints-1], c, c,
						this->parameters_t0[Npoints-1], this->parameters_t0[Npoints-1],
						parameters_fine[i], dparameters_fine[i]);
				bCheck = false;
			} else {
				for (int j=1; j<Npoints; ++j) { // Check at which point of the travelling wave we are
					if (( c >= positions[j-1] ) && ( c <= positions[j] )) {
						IF2D_Interpolation1D::cubicInterpolation(positions[j-1], positions[j], c,
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
		for(int i=Npoints-1; i>0; --i) this->parameters_t0[i] = this->parameters_t0[i-2];
		this->parameters_t0[1] = b;
		this->parameters_t0[0] = 0;
	}
};
}

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
			normalI[i]=std::abs(normalI[i])*invMagI;
			normalJ[i]=std::abs(normalJ[i])*invMagJ;
			normalK[i]=std::abs(normalK[i])*invMagK;
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

	bool isIntersectingWithAABB(const Real start[3],const Real end[3], const Real safe_distance = 0.0) const
	{
		//start and end are two diagonally opposed corners of grid block
		const Real AABB_w[3] = { //half block width + safe distance
				0.5*(end[0] - start[0]) + 2.0*safe_distance,
				0.5*(end[1] - start[1]) + 2.0*safe_distance,
				0.5*(end[2] - start[2]) + 2.0*safe_distance
		};
		const Real AABB_c[3] = { //block center
				start[0] + AABB_w[0] - safe_distance,
				start[1] + AABB_w[1] - safe_distance,
				start[2] + AABB_w[2] - safe_distance
		};
		assert(AABB_w[0]>0);
		assert(AABB_w[1]>0);
		assert(AABB_w[2]>0);
		bool intersects = true;
		Real r;
		{
			r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
			intersects &= ((c[0]-r <= AABB_c[0] + AABB_w[0]) && (c[0]+r >= AABB_c[0] - AABB_w[0]));

			r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
			intersects &= ((c[1]-r <= AABB_c[1] + AABB_w[1]) && (c[1]+r >= AABB_c[1] - AABB_w[1]));

			r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
			intersects &= ((c[2]-r <= AABB_c[2] + AABB_w[2]) && (c[2]+r >= AABB_c[2] - AABB_w[2]));
		}
		{
			//r = AABB_w[0]*normalI[0] + AABB_w[1]*normalI[1] + AABB_w[2]*normalI[2];
			r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
			intersects &= ((AABB_c[0]-r <= c[0] + w[0]) && (AABB_c[0]+r >= c[0] - w[0]));

			//r = AABB_w[0]*normalJ[0] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalJ[2];
			r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
			intersects &= ((AABB_c[1]-r <= c[1] + w[1]) && (AABB_c[1]+r >= c[1] - w[1]));

			//r = AABB_w[0]*normalK[0] + AABB_w[1]*normalK[1] + AABB_w[2]*normalK[2];
			r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
			intersects &= ((AABB_c[2]-r <= c[2] + w[2]) && (AABB_c[2]+r >= c[2] - w[2]));
		}
		return intersects;
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
			Real * const chi = &(defblock->chi[0][0][0]);
			Real * const udef = &(defblock->udef[0][0][0][0]);

			static const int n = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
			for(int i=0; i<n; i++) {
				chi[i]=-1;
				udef[3*i+0]=0;
				udef[3*i+1]=0;
				udef[3*i+2]=0;
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
																			  : (ss > cfish->iFishEnd   ? cfish->width[ cfish->iFishEnd] : cfish->width[ss]));
				const Real myHeight = (ss < cfish->iFishStart ? cfish->height[cfish->iFishStart]
																			  : (ss > cfish->iFishEnd   ? cfish->height[cfish->iFishEnd] : cfish->height[ss]));
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
									std::abs((Real)xp[c] - (ap[c]+0)),
									std::abs((Real)xp[c] - (ap[c]+1))
							};
							wghts[c][0] = 1.0 - t[0];
							wghts[c][1] = 1.0 - t[1];
						}

						const bool isInside = (std::abs(offsetW) < actualWidth) && (std::abs(offsetH) < cfish->height[ss]);
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
							if( (std::abs(defblock->chi[idx[2]][idx[1]][idx[0]] + 1) < 5*std::numeric_limits<Real>::epsilon()) && isInside)
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
		stencil = StencilInfo(-1,-1,-1, 2,2,2, true, 1, 5);
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}
#if 1
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
								.5*h+(Up[i]-fac1*sign(Up[i])*Up[i]*Up[i]));
				Hm[i] = (Um[i]> h) ? h : (
						(Um[i]<-h) ? 0 :
								.5*h+(Um[i]-fac1*sign(Um[i])*Um[i]*Um[i]));
			}
			const Real gradH[3] = {.5*(Hp[0]-Hm[0]), .5*(Hp[1]-Hm[1]), .5*(Hp[2]-Hm[2])};
			Real gradUU[3], gradUI[3], gradUH[3];
			for (int i=0; i<3; i++) {
				gradUU[i] = gradU[i]*gradU[i] + gradU[i+3]*gradU[i+3] + gradU[i+6]*gradU[i+6];
				gradUI[i] = gradU[i]*gradI[i] + gradU[i+3]*gradI[i+3] + gradU[i+6]*gradI[i+6];
				gradUH[i] = gradU[i]*gradH[i];
			}
			for (int i=0; i<3; i++)  gradUU[i] = max(gradUU[i], eps);
			const Real FDD = h*(gradUH[0] + gradUH[1] + gradUH[2])/gradUU[0]; // == delta * h^3
			const Real FDH = 1/3. * (gradUI[0]/gradUU[0]+gradUI[1]/gradUU[1]+gradUI[2]/gradUU[2]);


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

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU >= +2*h || lab(ix,iy,iz).tmpU <= -2*h) {
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

			// gradU
			const Real gradUX = inv2h * (distPx - distMx);
			const Real gradUY = inv2h * (distPy - distMy);
			const Real gradUZ = inv2h * (distPz - distMz);
			const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ;

			// gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
			const Real gradIX = inv2h * (IplusX - IminuX);
			const Real gradIY = inv2h * (IplusY - IminuY);
			const Real gradIZ = inv2h * (IplusZ - IminuZ);
			const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;

			const Real gradHX = inv2h * (HplusX - HminuX);
			const Real gradHY = inv2h * (HplusY - HminuY);
			const Real gradHZ = inv2h * (HplusZ - HminuZ);
			const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;

			const Real Delta = std::abs(gradUSq) < eps ? numD : numD/gradUSq;
			const Real H     = std::abs(gradUSq) < eps ? numH : numH/gradUSq;

			if (Delta>1e-6) {
				const Real dchidx = -Delta*gradUX;
				const Real dchidy = -Delta*gradUY;
				const Real dchidz = -Delta*gradUZ;
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, Delta);
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

	CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext);
	void _correctTrajectory(const Real dtheta, const Real time, const Real dt) override
			void execute(const Real time, const Real l_tnext, const vector<Real>& input) override;
	~CurvatureDefinedFishData();
	void computeMidline(const Real time);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
