//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_FishLibrary__
#define __IncompressibleFluids3D__IF3D_FishLibrary__
//#define BBURST
#include <cmath>
#include <array>

#include "Definitions.h"
#include "IF2D_Interpolation1D.h"
#include "GenericOperator.h"
#include "IF2D_Frenet.h"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>

const int NPPEXT = 2; //was 3, BUT now i do towers only?
const int TGTPPB = 2; //was 2 i think
const int NEXTDX = 2; //was 4
#define __BSPLINE

namespace Schedulers
{
template<int Npoints>
struct ParameterScheduler
{
  std::array<double, Npoints> parameters_t0; // parameters at t0
  std::array<double, Npoints> parameters_t1; // parameters at t1
  std::array<double, Npoints> dparameters_t0; // derivative at t0
  double t0, t1; // t0 and t1

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
    parameters_t0 = std::array<double, Npoints>();
    parameters_t1 = std::array<double, Npoints>();
    dparameters_t0 = std::array<double, Npoints>();
  }

  void transition(const double t, const double tstart, const double tend,
      const std::array<double, Npoints> parameters_tend,
      const bool UseCurrentDerivative = false)
  {
    if(t<tstart or t>tend) return; // this transition is out of scope
    //if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

    // we transition from whatever state we are in to a new state
    // the start point is where we are now: lets find out
    std::array<double, Npoints> parameters;
    std::array<double, Npoints> dparameters;
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
    dparameters_t0 = UseCurrentDerivative ? dparameters : std::array<double, Npoints>();
  }

  void transition(const double t, const double tstart, const double tend,
      const std::array<double, Npoints> parameters_tstart,
      const std::array<double, Npoints> parameters_tend)
  {
    if(t<tstart or t>tend) return; // this transition is out of scope
    if(tstart<t0) return; // this transition is not relevant: we are doing a next one already

    // fill my members
    t0 = tstart;
    t1 = tend;
    parameters_t0 = parameters_tstart;
    parameters_t1 = parameters_tend;
  }

  void gimmeValues(const double t, std::array<double, Npoints>& parameters, std::array<double, Npoints>& dparameters)
  {
    // look at the different cases
    if(t<t0 or t0<0) { // no transition, we are in state 0
      parameters = parameters_t0;
      dparameters = std::array<double, Npoints>();
    } else if(t>t1) { // no transition, we are in state 1
      parameters = parameters_t1;
      dparameters = std::array<double, Npoints>();
    } else { // we are within transition: interpolate
      for(int i=0;i<Npoints;++i)
        IF2D_Interpolation1D::cubicInterpolation(t0,t1,t,parameters_t0[i],parameters_t1[i],dparameters_t0[i],0.0,parameters[i],dparameters[i]);
    }
  }

  void gimmeValues(const double t, std::array<double, Npoints>& parameters)
  {
    std::array<double, Npoints> dparameters_whocares; // no derivative info
    return gimmeValues(t,parameters,dparameters_whocares);
  }
};

struct ParameterSchedulerScalar : ParameterScheduler<1>
{
  void transition(const double t, const double tstart, const double tend, const double parameter_tend, const bool UseCurrentDerivative = false)
  {
    const std::array<double, 1> myParameter = {parameter_tend};
    return ParameterScheduler<1>::transition(t,tstart,tend,myParameter,UseCurrentDerivative);
  }

  void gimmeValues(const double t, double & parameter, double & dparameter)
  {
    std::array<double, 1> myParameter, mydParameter;
    ParameterScheduler<1>::gimmeValues(t, myParameter, mydParameter);
    parameter = myParameter[0];
    dparameter = mydParameter[0];
  }

  void gimmeValues(const double t, double & parameter)
  {
    std::array<double, 1> myParameter;
    ParameterScheduler<1>::gimmeValues(t, myParameter);
    parameter = myParameter[0];
  }
};

template<int Npoints>
struct ParameterSchedulerVector : ParameterScheduler<Npoints>
{
  void gimmeValues(const double t, const std::array<double, Npoints> & positions, const int Nfine,
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
            0.0,           parameters_fine[i],   dparameters_fine[i]);
    }
  }

  void gimmeValues(const double t, std::array<double, Npoints>& parameters)
  {
    ParameterScheduler<Npoints>::gimmeValues(t, parameters);
  }

  void gimmeValues(const double t, std::array<double, Npoints> & parameters, std::array<double, Npoints> & dparameters)
  {
    ParameterScheduler<Npoints>::gimmeValues(t, parameters, dparameters);
  }
};

template<int Npoints>
struct ParameterSchedulerLearnWave : ParameterScheduler<Npoints>
{
  void gimmeValues(const double t, const double Twave, const double Length,
    const std::array<double, Npoints> & positions, const int Nfine,
    const Real* const positions_fine, Real* const parameters_fine, Real* const dparameters_fine)
  {
    const double _1oL = 1./Length;
    const double _1oT = 1./Twave;
    // the fish goes through (as function of t and s) a wave function that describes the curvature
    for(int i=0;i<Nfine;++i) {
      const double c = positions_fine[i]*_1oL - (t - this->t0)*_1oT; //traveling wave coord
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

  void Turn(const double b, const double t_turn) // each decision adds a node at the beginning of the wave (left, right, straight) and pops last node
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
  const double length, Tperiod, phaseShift;
  double l_Tp, timeshift, time0;
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
  const int iFishStart, iFishEnd;
  Real * const rXold;
  Real * const rYold;
  Real oldTime = 0.0;
  bool firstStep = true;

  double linMom[2], vol, J, angMom; // for diagnostics
  // start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
  Schedulers::ParameterSchedulerVector<6> curvScheduler;
  Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
  Schedulers::ParameterSchedulerVector<6> adjustScheduler;
  FishSkin * upperSkin, * lowerSkin;
  //Real finSize = 1.1e-1, waveLength = 1.0;
  Real finSize = 0.0e-1, waveLength = 1.0; // For curvalicious fin

 protected:
  double Rmatrix2D[2][2];
  double Rmatrix3D[3][3];

  inline void _rotate2D(Real &x, Real &y) const
  {
    const double p[2] = {x,y};
    x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
    y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
  }

  inline void _translateAndRotate2D(const Real pos[2], Real &x, Real &y) const
  {
    const double p[2] = {
        x-pos[0],
        y-pos[1]
    };
    // rotate
    x = Rmatrix2D[0][0]*p[0] + Rmatrix2D[0][1]*p[1];
    y = Rmatrix2D[1][0]*p[0] + Rmatrix2D[1][1]*p[1];
  }

  inline double _d_ds(const int idx, const Real* const vals, const int maxidx) const
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

  inline double _integrationFac1(const int idx) const
  {
    return double(width[idx])*height[idx];
  }

  inline double _integrationFac2(const int idx) const
  {
    const double dnorXi = _d_ds(idx, norX, Nm);
    const double dnorYi = _d_ds(idx, norY, Nm);
    return 0.25*std::pow((double)width[idx],3)*height[idx]*(dnorXi*norY[idx] - dnorYi*norX[idx]);
  }

  inline double _integrationFac3(const int idx) const
  {
    // const double drXi = _d_ds(idx, rX, Nm);
    // const double drYi = _d_ds(idx, rY, Nm);
    // return 0.25*std::pow(width[idx],3)*height[idx]*(drXi*norY[idx] - drYi*norX[idx]);
    return 0.25*std::pow((double)width[idx],3)*height[idx];
  }

  void _prepareRotation2D(const double angle)
  {
    Rmatrix2D[0][0] = Rmatrix2D[1][1] = std::cos(angle);
    Rmatrix2D[0][1] = -std::sin(angle);
    Rmatrix2D[1][0] = -Rmatrix2D[0][1];
  }

  /*
  function inputs: xc, yc are n sized arrays which contain the control points of the cubic b spline
  function outputs onto res: assumed to be either the width or the height
  */
  void integrateBSpline(Real* const res, const double* const xc,
                                         const double* const yc, const int n);

  void _computeMidlineNormals();
  virtual void computeWidthsHeights() = 0;

 public:
  FishMidlineData(const int _Nm, const double len, const double Tp, const double phase, const double dx_ext):
   Nm(_Nm),length(len),Tperiod(Tp),phaseShift(phase),l_Tp(Tperiod),timeshift(0),
   time0(0), rS(_alloc(_Nm)), rX(_alloc(_Nm)), rY(_alloc(_Nm)), vX(_alloc(_Nm)),
   vY(_alloc(_Nm)), norX(_alloc(_Nm)), norY(_alloc(_Nm)), vNorX(_alloc(_Nm)),
   vNorY(_alloc(_Nm)), width(_alloc(_Nm)), height(_alloc(_Nm)),
   iFishStart(NEXTDX*NPPEXT), iFishEnd(_Nm-1-NEXTDX*NPPEXT),
   rXold(_alloc(_Nm)), rYold(_alloc(_Nm))
  {
    std::fill(rXold, rXold+Nm, 0.0);
    std::fill(rYold, rYold+Nm, 0.0);

    // extension_info contains number of extension points and extension dx
    const int Nextension = NEXTDX*NPPEXT; // up to 3dx on each side (to get proper interpolation up to 2dx)
    const int Next = Nextension; // number of points per extension
    const int Nint = Nm -2*Next; // number of interior points

    // extension head
    for(int i=0;i<Next;++i) rS[i] = 0 - (Next-i)*dx_ext;
    // interior points
    for(int i=0;i<Nint;++i)
      rS[i+Next] = length * 0.5 * (1.0 - std::cos(i * M_PI/((Real)Nint-1)));
      // cosine: more points near head and tail
    // rS[i] = i*length/((Real)Nint-1); // linear: equally distributed points
    // extension tail
    for(int i=0;i<Next;++i) rS[i+Nint+Next] = length +(i+1)*dx_ext;
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

  Real integrateLinearMomentum(double CoM[2], double vCoM[2]);

  void integrateAngularMomentum(double& angVel);

  void changeToCoMFrameLinear(const double CoM_internal[2], const double vCoM_internal[2]);

  void changeToCoMFrameAngular(const double theta_internal, const double angvel_internal);

  void computeSurface();

  void computeSkinNormals(const double theta_comp, const double CoM_comp[3]);

  void surfaceToCOMFrame(const double theta_internal, const double CoM_internal[2]);

  void surfaceToComputationalFrame(const double theta_comp, const double CoM_interpolated[3]);

  virtual void computeMidline(const double time) = 0;

  virtual void _correctAmplitude(const double dAmp, const double vAmp, const double time, const double dt) {}
  virtual void _correctTrajectory(const double dtheta, const double vtheta, const double time, const double dt) {}
  virtual void execute(const double time, const double l_tnext, const vector<double>& input) {}
};

struct VolumeSegment_OBB
{
  std::pair<int, int> s_range;
  Real normalI[3] = {1,0,0}; // should be normalized and >=0
  Real normalJ[3] = {0,1,0};
  Real normalK[3] = {0,0,1};
  Real w[3]; // halfwidth
  Real c[3]; // center

  VolumeSegment_OBB() { }

  void prepare(std::pair<int, int> _s_range, const Real bbox[3][2]);

  void normalizeNormals();

  void changeToComputationalFrame(const double position[3], const double quaternion[4]);

  bool isIntersectingWithAABB(const Real start[3],const Real end[3], const Real safe_distance = 0.0) const;
};

struct PutFishOnBlocks
{
  const FishMidlineData * cfish;
  const double position[3];
  const double quaternion[4];
  const double Rmatrix3D[3][3];

  PutFishOnBlocks(const FishMidlineData* const cfish, const double p[3], const double q[4]): cfish(cfish), position{p[0],p[1],p[2]},
  quaternion{q[0],q[1],q[2],q[3]},
  Rmatrix3D{
  {1-2*(q[2]*q[2]+q[3]*q[3]), 2*(q[1]*q[2]-q[3]*q[0]), 2*(q[1]*q[3]+q[2]*q[0])},
  {2*(q[1]*q[2]+q[3]*q[0]), 1-2*(q[1]*q[1]+q[3]*q[3]), 2*(q[2]*q[3]-q[1]*q[0])},
  {2*(q[1]*q[3]-q[2]*q[0]), 2*(q[2]*q[3]+q[1]*q[0]), 1-2*(q[1]*q[1]+q[2]*q[2])}
  } { }

  inline int find_closest_dist(const int s, const int dir, const Real x[3], Real & oldDistSq) const
  {
    if((s+dir)<cfish->iFishStart or (s+dir)>cfish->iFishEnd) return s;

    const Real newDistSq = (x[0]-cfish->rX[s+dir])*(x[0]-cfish->rX[s+dir])
      + (x[1]-cfish->rY[s+dir])*(x[1]-cfish->rY[s+dir]) + (x[2])*(x[2]);
    if(oldDistSq<=newDistSq) return s;
    else {
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
    const T p[3] = { // translate back to CoM
        x[0]-position[0],
        x[1]-position[1],
        x[2]-position[2]
    };
    // rotate back around CoM
    x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[1][0]*p[1] + Rmatrix3D[2][0]*p[2];
    x[1]=Rmatrix3D[0][1]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[2][1]*p[2];
    x[2]=Rmatrix3D[0][2]*p[0] + Rmatrix3D[1][2]*p[1] + Rmatrix3D[2][2]*p[2];
  }

  Real getSmallerDistToMidline(const int start_s, const Real x[3], int & final_s) const;

  void operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const;
};

struct PutFishOnBlocks_Finalize : public GenericLabOperator
{
  Real t;
  const int stencil_start[3] = {-1, -1, -1}, stencil_end[3] = {2, 2, 2};
  const Real eps = std::numeric_limits<Real>::epsilon();

  PutFishOnBlocks_Finalize() : t(0)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 5);
  }

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& b, ObstacleBlock*const o)
  {
    const Real h = info.h_gridpoint;
    const Real inv2h = .5/h;
    const Real fac1 = 0.5*h*h;
    const Real fac2 = h*h*h;

    for(int iz=0; iz<FluidBlock::sizeZ; iz++)
    for(int iy=0; iy<FluidBlock::sizeY; iy++)
    for(int ix=0; ix<FluidBlock::sizeX; ix++) {
      Real p[3];
      info.pos(p, ix,iy,iz);
      if (lab(ix,iy,iz).tmpU > +2*h || lab(ix,iy,iz).tmpU < -2*h)
      {
        const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
        b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
        o->write(ix,iy,iz,H,0,0,0,0,0);
        o->CoM_x += p[0]*H;
        o->CoM_y += p[1]*H;
        o->CoM_z += p[2]*H;
        o->mass += H;
        continue;
      }

      const Real distPx = lab(ix+1,iy,iz).tmpU, distMx = lab(ix-1,iy,iz).tmpU;
      const Real distPy = lab(ix,iy+1,iz).tmpU, distMy = lab(ix,iy-1,iz).tmpU;
      const Real distPz = lab(ix,iy,iz+1).tmpU, distMz = lab(ix,iy,iz-1).tmpU;
      // gradU
      const Real gradUX = inv2h*(distPx - distMx);
      const Real gradUY = inv2h*(distPy - distMy);
      const Real gradUZ = inv2h*(distPz - distMz);
      const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ + eps;

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

      o->write(ix, iy, iz, H, Delta, gradUX, gradUY, gradUZ, fac1);
      b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
      o->CoM_x += p[0]*H;
      o->CoM_y += p[1]*H;
      o->CoM_z += p[2]*H;
      o->mass += H;
    }
  }
};

#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
