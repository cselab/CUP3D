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
#include "GenericOperator.h"
#include "IF2D_Frenet.h"
#include "IF3D_Schedulers.h"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>

#define __BSPLINE

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
  const double length, Tperiod, phaseShift, h;
  const Real waveLength = 1;

  // Midline is discretized by more points in first fraction and last fraction:
  const double fracRefined = 0.1, fracMid = 1 - 2*fracRefined;
  const double dSmid_tgt = h / std::sqrt(3);
  const double dSrefine_tgt = 0.125 * h;

  const int Nmid = (int)std::floor(length * fracMid / dSmid_tgt / 8) * 8;
  const double dSmid = length * fracMid / Nmid;

  const int Nend = (int)std::ceil( // here we ceil to be safer
    fracRefined * length * 2 / (dSmid + dSrefine_tgt)  / 8) * 8;
  const double dSref = fracRefined * length * 2 / Nend - dSmid;

  const int Nm = Nmid + 2 * Nend + 1; // plus 1 because we contain 0 and L

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
  double * const forceX;
  double * const forceY;
  double * const torque;
  Real oldTime = 0.0;
  // quantities needed to correctly control the speed of the midline maneuvers
  double l_Tp = Tperiod, timeshift = 0, time0 = 0;
  bool firstStep = true;

  double linMom[2], vol, J, angMom; // for diagnostics
  // start and end indices in the arrays where the fish starts and ends (to ignore the extensions when interpolating the shapes)
  Schedulers::ParameterSchedulerVector<6> curvScheduler;
  Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
  Schedulers::ParameterSchedulerVector<6> adjustScheduler;
  FishSkin * upperSkin, * lowerSkin;

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

  template<typename T>
  void _dealloc(T * ptr)
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

  void _computeMidlineNormals();

 public:
  FishMidlineData(double L, double Tp, double phi, double _h):
   length(L), Tperiod(Tp), phaseShift(phi), h(_h), rS(_alloc(Nm)),
   rX(_alloc(Nm)), rY(_alloc(Nm)), vX(_alloc(Nm)), vY(_alloc(Nm)),
   norX(_alloc(Nm)), norY(_alloc(Nm)), vNorX(_alloc(Nm)), vNorY(_alloc(Nm)),
   width(_alloc(Nm)), height(_alloc(Nm)), forceX(new double[Nm]),
   forceY(new double[Nm]), torque(new double[Nm]),
   upperSkin(new FishSkin(Nm)), lowerSkin(new FishSkin(Nm))
  {
    std::fill(forceX, forceX+Nm, 0.0); // these are initialized to 0 because
    std::fill(forceY, forceY+Nm, 0.0); // force compute is at end of time step
    std::fill(torque, torque+Nm, 0.0);
    // extension head
    rS[0] = 0;
    int k = 0;
    for(int i=0; i<Nend; ++i, k++)
      rS[k+1] = rS[k] + dSref +(dSmid-dSref) *         i /((Real)Nend-1.);

    // interior points
    for(int i=0; i<Nmid; ++i, k++)
      rS[k+1] = rS[k] + dSmid;

    // extension tail
    for(int i=0; i<Nend; ++i, k++)
      rS[k+1] = rS[k] + dSref +(dSmid-dSref) * (Nend-i-1)/((Real)Nend-1.);

    assert(k+1==Nm);
    //cout << "Discrepancy of midline length: " << std::fabs(rS[k]-L) << endl;
    rS[k] = std::min(rS[k], (Real)L);
  }

  virtual ~FishMidlineData()
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
    _dealloc(forceX);
    _dealloc(forceY);
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

  virtual void computeMidline(const double time, const double dt) = 0;

  virtual void _correctAmplitude(const double dAmp, const double vAmp, const double time, const double dt) {}
  virtual void _correctTrajectory(const double dtheta, const double vtheta, const double time, const double dt) {}
  virtual void execute(const double time, const double l_tnext, const vector<double>& input) {}
  void writeMidline2File(const int step_id, std::string filename);
};

struct VolumeSegment_OBB
{
  Real safe_distance = 0;
  std::pair<int, int> s_range;
  Real normalI[3] = {1,0,0}; // should be normalized and >=0
  Real normalJ[3] = {0,1,0};
  Real normalK[3] = {0,0,1};
  Real w[3]={0,0,0}, c[3]={0,0,0}; // halfwidth & center
  Real objBoxLabFr[3][2] = {{0,0}, {0,0}, {0,0}};
  Real objBoxObjFr[3][2] = {{0,0}, {0,0}, {0,0}};

  VolumeSegment_OBB() { }

  void prepare(std::pair<int, int> _s_range, const Real bbox[3][2], const Real safe_dist);

  void normalizeNormals();

  void changeToComputationalFrame(const double position[3], const double quaternion[4]);

  bool isIntersectingWithAABB(const Real start[3],const Real end[3]) const;
};

struct PutFishOnBlocks
{
  const FishMidlineData * cfish;
  const double position[3];
  const double quaternion[4];
  const double Rmatrix3D[3][3];

  PutFishOnBlocks(const FishMidlineData* const _cfish, const double p[3], const double q[4]): cfish(_cfish), position{p[0],p[1],p[2]},
  quaternion{q[0],q[1],q[2],q[3]},
  Rmatrix3D{
  {1-2*(q[2]*q[2]+q[3]*q[3]), 2*(q[1]*q[2]-q[3]*q[0]), 2*(q[1]*q[3]+q[2]*q[0])},
  {2*(q[1]*q[2]+q[3]*q[0]), 1-2*(q[1]*q[1]+q[3]*q[3]), 2*(q[2]*q[3]-q[1]*q[0])},
  {2*(q[1]*q[3]-q[2]*q[0]), 2*(q[2]*q[3]+q[1]*q[0]), 1-2*(q[1]*q[1]+q[2]*q[2])}
  } { }

  static inline Real eulerDistSq3D(const Real a[3], const Real b[3]) {
    return std::pow(a[0]-b[0],2) +std::pow(a[1]-b[1],2) +std::pow(a[2]-b[2],2);
  }
  static inline Real eulerDistSq2D(const Real a[3], const Real b[3]) {
    return std::pow(a[0]-b[0],2) +std::pow(a[1]-b[1],2);
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
        x[0]-(T)position[0],
        x[1]-(T)position[1],
        x[2]-(T)position[2]
    };
    // rotate back around CoM
    x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[1][0]*p[1] + Rmatrix3D[2][0]*p[2];
    x[1]=Rmatrix3D[0][1]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[2][1]*p[2];
    x[2]=Rmatrix3D[0][2]*p[0] + Rmatrix3D[1][2]*p[1] + Rmatrix3D[2][2]*p[2];
  }

  void operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const;
  virtual void constructSurface(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const oblck, const std::vector<VolumeSegment_OBB>& vSeg) const;
  virtual void constructInternl(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const oblck, const std::vector<VolumeSegment_OBB>& vSeg) const;
  virtual void signedDistanceSqrt(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const oblck, const std::vector<VolumeSegment_OBB>& vSeg) const;
};

struct PutNacaOnBlocks: public PutFishOnBlocks
{
  PutNacaOnBlocks(const FishMidlineData* const _cfish, const double p[3], const double q[4]): PutFishOnBlocks(_cfish, p, q) { }

  Real getSmallerDistToMidLPlanar(const int start_s, const Real x[3], int & final_s) const;

  void constructSurface(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegm) const override;
  void constructInternl(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegm) const override;
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
    static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
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
      const Real HplusX = std::fabs(distPx)<EPS ? 0.5 : (distPx < 0 ? 0 : 1);
      const Real HminuX = std::fabs(distMx)<EPS ? 0.5 : (distMx < 0 ? 0 : 1);
      const Real HplusY = std::fabs(distPy)<EPS ? 0.5 : (distPy < 0 ? 0 : 1);
      const Real HminuY = std::fabs(distMy)<EPS ? 0.5 : (distMy < 0 ? 0 : 1);
      const Real HplusZ = std::fabs(distPz)<EPS ? 0.5 : (distPz < 0 ? 0 : 1);
      const Real HminuZ = std::fabs(distMz)<EPS ? 0.5 : (distMz < 0 ? 0 : 1);

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

namespace MidlineShapes
{
  /*
  function inputs: xc, yc are n sized arrays which contain the control points of the cubic b spline
  function outputs onto res: assumed to be either the width or the height
  */
  void integrateBSpline(const double*const xc, const double*const yc,
  const int n, const double length, Real*const rS,Real*const res,const int Nm);

  void naca_width(const double t_ratio, const double L, Real*const rS,
    Real*const res, const int Nm);
  void stefan_width(const double L, Real*const rS, Real*const res, const int Nm);
  void stefan_height(const double L, Real*const rS, Real*const res, const int Nm);
  void danio_width(const double L, Real*const rS, Real*const res, const int Nm);
  void danio_height(const double L, Real*const rS, Real*const res, const int Nm);

  void computeWidthsHeights(const string heightName, const string widthName,
  const double L, Real*const rS, Real*const height, Real*const width,
  const int nM, const int mpirank);
}
#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
