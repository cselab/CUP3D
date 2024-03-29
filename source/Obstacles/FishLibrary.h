//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "../Definitions.h"
#include "../ObstacleBlock.h"
#include "extra/Frenet.h"
#include "extra/Schedulers.h"

#include <cmath>

CubismUP_3D_NAMESPACE_BEGIN

class FishMidlineData
{
 public:
  const Real length;     //midline length
  const Real Tperiod;    //tail-beat period
  const Real phaseShift; //phase shift for tail-beat
  const Real h;          //grid spacing of grid where midline will be created
  const Real waveLength = 1; // midline (or its curvature) parametrized as wave with this wavelength
  const Real amplitudeFactor;// midline curvature amplitude

  //Midline discretization parameters
  //Midline is discretized by more points in first fraction and last fraction:
  const Real fracRefined = 0.1;
  const Real fracMid = 1 - 2*fracRefined;
  const Real dSmid_tgt = h / std::sqrt(3);
  const Real dSrefine_tgt = 0.125 * h;
  const int Nmid = (int)std::ceil(length * fracMid / dSmid_tgt / 8) * 8;
  const Real dSmid = length * fracMid / Nmid;
  const int Nend = (int)std::ceil(fracRefined * length * 2 / (dSmid + dSrefine_tgt) / 4) * 4;
  const Real dSref = fracRefined * length * 2 / Nend - dSmid;
  const int Nm = Nmid + 2 * Nend + 1; // plus 1 because we contain 0 and L

  Real * const rS;   // arclength discretization points

  Real * const rX;   // X-coordinate of midline discretization points
  Real * const rY;   // Y-coordinate of midline discretization points
  Real * const rZ;   // Z-coordinate of midline discretization points
  Real * const vX;   // midline discretization velocities (=drX/dt)
  Real * const vY;   // midline discretization velocities (=drY/dt)
  Real * const vZ;   // midline discretization velocities (=drZ/dt)

  Real * const norX; // normal vector to the midline discretization points (X component)
  Real * const norY; // normal vector to the midline discretization points (Y component)
  Real * const norZ; // normal vector to the midline discretization points (Z component)
  Real * const vNorX;// time derivative of normal vector (=dnorX/dt)
  Real * const vNorY;// time derivative of normal vector (=dnorY/dt)
  Real * const vNorZ;// time derivative of normal vector (=dnorZ/dt)

  Real * const binX; // binormal vector to the midline discretization points (X component)
  Real * const binY; // binormal vector to the midline discretization points (Y component)
  Real * const binZ; // binormal vector to the midline discretization points (Z component)
  Real * const vBinX;// time derivative of binormal vector (=dbinX/dt)
  Real * const vBinY;// time derivative of binormal vector (=dbinY/dt)
  Real * const vBinZ;// time derivative of binormal vector (=dbinZ/dt)

  //fish cross-section is an ellipse with axis 2*width(s) and 2*height(s)
  //the surface is parametrized by two parameters: 
  // s    : coordinate along midline (discretized by rS), in (0,L)
  // theta: angle in a cross-section, in (0,2pi)
  // The parametrization is: x(s,theta) = rS(s) + nor(s)*width(s)*cos(theta) + bin(s)*height(s)*sin(theta)
  Real * const width;
  Real * const height;

  std::array<Real, 9> sensorLocation; //Shear stress sensor locations (for RL)

  // Midline has an orientation in space which is defined from the following quaternion.
  // When midline is defined, we change the frame of reference so that its origin is at the midline
  // center of mass and its linear and angular momentums are zero.
  Real quaternion_internal[4]={1,0,0,0};
  Real angvel_internal[3]={0,0,0};

 protected:

  //computes derivative d/ds of given quantity
  inline Real _d_ds(const int idx, const Real* const vals, const int maxidx) const
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

 public:
  FishMidlineData(Real L, Real Tp, Real phi, Real _h, Real _ampFac=1):
   length(L), Tperiod(Tp), phaseShift(phi), h(_h), amplitudeFactor(_ampFac),
   rS    (_alloc(Nm)),
   rX    (_alloc(Nm)), rY   (_alloc(Nm)), rZ   (_alloc(Nm)),
   vX    (_alloc(Nm)), vY   (_alloc(Nm)), vZ   (_alloc(Nm)),
   norX  (_alloc(Nm)), norY (_alloc(Nm)), norZ (_alloc(Nm)),
   vNorX (_alloc(Nm)), vNorY(_alloc(Nm)), vNorZ(_alloc(Nm)),
   binX  (_alloc(Nm)), binY (_alloc(Nm)), binZ (_alloc(Nm)),
   vBinX (_alloc(Nm)), vBinY(_alloc(Nm)), vBinZ(_alloc(Nm)),
   width (_alloc(Nm)), height(_alloc(Nm))
  {
    // Define points along midline
    rS[0] = 0;
    int k = 0;
    for(int i=0; i<Nend; ++i, k++) //fish head
      rS[k+1] = rS[k] + dSref +(dSmid-dSref) *         i /((Real)Nend-1.);
    for(int i=0; i<Nmid; ++i, k++) //interior points
      rS[k+1] = rS[k] + dSmid;
    for(int i=0; i<Nend; ++i, k++) //fish tail
      rS[k+1] = rS[k] + dSref +(dSmid-dSref) * (Nend-i-1)/((Real)Nend-1.);
    rS[k] = std::min(rS[k], (Real)L);
    assert(k+1==Nm);
  }

  virtual ~FishMidlineData()
  {
    _dealloc(rS);
    _dealloc(   rX); _dealloc(   rY); _dealloc(   rZ);
    _dealloc(   vX); _dealloc(   vY); _dealloc(   vZ);
    _dealloc( norX); _dealloc( norY); _dealloc( norZ);
    _dealloc(vNorX); _dealloc(vNorY); _dealloc(vNorZ);
    _dealloc( binX); _dealloc( binY); _dealloc( binZ);
    _dealloc(vBinX); _dealloc(vBinY); _dealloc(vBinZ);
    _dealloc(width);
    _dealloc(height);
  }

  void writeMidline2File(const int step_id, std::string filename)
  {
    char buf[500];
    sprintf(buf, "%s_midline_%07d.txt", filename.c_str(), step_id);
    FILE * f = fopen(buf, "a");
    fprintf(f, "s x y z vX vY vZ\n");
    for (int i=0; i<Nm; i++)
      fprintf(f, "%g %g %g %g %g %g %g\n", rS[i],rX[i],rY[i],rZ[i],vX[i],vY[i],vZ[i]);
    fclose(f);
  }

  void writeMidline2File(const int step_id, std::string filename, const double q[4], const double position[3])
  {
    const double Rmatrix3D[3][3] = {
    {1-2*(q[2]*q[2]+q[3]*q[3]), 2*(q[1]*q[2]-q[3]*q[0]), 2*(q[1]*q[3]+q[2]*q[0])},
    {2*(q[1]*q[2]+q[3]*q[0]), 1-2*(q[1]*q[1]+q[3]*q[3]), 2*(q[2]*q[3]-q[1]*q[0])},
    {2*(q[1]*q[3]-q[2]*q[0]), 2*(q[2]*q[3]+q[1]*q[0]), 1-2*(q[1]*q[1]+q[2]*q[2])}};

    char buf[500];
    sprintf(buf, "%s_midline_%07d.txt", filename.c_str(), step_id);
    FILE * f = fopen(buf, "a");
    fprintf(f, "x y z s\n");
    for (int i=0; i<Nm; i++)
    {
      double x[3];
      x[0]=position[0] + Rmatrix3D[0][0]*rX[i] + Rmatrix3D[0][1]*rY[i] + Rmatrix3D[0][2]*rZ[i];
      x[1]=position[1] + Rmatrix3D[1][0]*rX[i] + Rmatrix3D[1][1]*rY[i] + Rmatrix3D[1][2]*rZ[i];
      x[2]=position[2] + Rmatrix3D[2][0]*rX[i] + Rmatrix3D[2][1]*rY[i] + Rmatrix3D[2][2]*rZ[i];
      fprintf(f, "%g %g %g %g \n", x[0],x[1],x[2],rS[i]);
    }
    fclose(f);
  }



  void SurfaceNormal(const int idx, const Real theta, Real & nx, Real & ny, Real &nz)
  {
    //Surface is given from:
    //x(s,theta) = rS(s) + nor(s)*width(s)*cos(theta) + bin(s)*height(s)*sin(theta)
    //Here we compute n := dx/dtheta x dx/ds and we normalize the result, 
    //                                       to get the unit normal vector to the fish surface

    //if s = 0 then width = 0, height = 0 and the function is not differentiable.
    //in this case, we DEFINE the "normal" vector as the tangent vector dr/ds 
    //(i.e. the fish is looking forward)  
    if (idx == 0)
    {
      nx = (rX[idx+1]-rX[idx])/(rS[idx+1]-rS[idx]);
      ny = (rY[idx+1]-rY[idx])/(rS[idx+1]-rS[idx]);
      nz = (rZ[idx+1]-rZ[idx])/(rS[idx+1]-rS[idx]);
    }
    else
    {
      const Real costheta = cos(theta);
      const Real sintheta = sin(theta);

      Real drXds   = _d_ds(idx,  rX,Nm);
      Real drYds   = _d_ds(idx,  rY,Nm);
      Real drZds   = _d_ds(idx,  rZ,Nm);
      Real dnorXds = _d_ds(idx,norX,Nm);
      Real dnorYds = _d_ds(idx,norY,Nm);
      Real dnorZds = _d_ds(idx,norZ,Nm);
      Real dbinXds = _d_ds(idx,binX,Nm);
      Real dbinYds = _d_ds(idx,binY,Nm);
      Real dbinZds = _d_ds(idx,binZ,Nm);
      Real dwds    = _d_ds(idx, width,Nm);
      Real dhds    = _d_ds(idx,height,Nm);

      Real dxds = drXds + costheta * (dwds*norX[idx] + dnorXds*width[idx]) + sintheta * (dhds*binX[idx]+dbinXds*height[idx]);
      Real dyds = drYds + costheta * (dwds*norY[idx] + dnorYds*width[idx]) + sintheta * (dhds*binY[idx]+dbinYds*height[idx]);
      Real dzds = drZds + costheta * (dwds*norZ[idx] + dnorZds*width[idx]) + sintheta * (dhds*binZ[idx]+dbinZds*height[idx]);

      Real dxdtheta = -sintheta*norX[idx]*width[idx] + costheta*binX[idx]*height[idx];
      Real dydtheta = -sintheta*norY[idx]*width[idx] + costheta*binY[idx]*height[idx];
      Real dzdtheta = -sintheta*norZ[idx]*width[idx] + costheta*binZ[idx]*height[idx];

      nx = dydtheta * dzds - dzdtheta * dyds;
      ny = dzdtheta * dxds - dxdtheta * dzds;
      nz = dxdtheta * dyds - dydtheta * dxds;
    }
    const Real norm = 1.0/(sqrt(nx*nx+ny*ny+nz*nz)+std::numeric_limits<Real>::epsilon());
    nx *= norm;
    ny *= norm;
    nz *= norm;
  }

  void integrateLinearMomentum();

  void integrateAngularMomentum(const Real dt);

  // Derived class should provide the following function, which is responsible for defining 
  // rX,rY,rZ,norX,norY,norZ,binX,binY,biZ,vX,vY,vZ,vNorX,vNorY,vNorZ,vBinX,vBinY,vBinZ 
  virtual void computeMidline(const Real time, const Real dt) = 0;

  // used in RL
  virtual void execute(const Real time, const Real l_tnext, const std::vector<Real>& input) {}
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

  void changeToComputationalFrame(const Real position[3], const Real quaternion[4]);

  bool isIntersectingWithAABB(const Real start[3],const Real end[3]) const;
};

struct PutFishOnBlocks
{
  FishMidlineData * cfish;
  const Real position[3];
  const Real quaternion[4];
  const Real Rmatrix3D[3][3];

  PutFishOnBlocks(FishMidlineData* _cfish, const Real p[3], const Real q[4]): cfish(_cfish), position{p[0],p[1],p[2]},
  quaternion{q[0],q[1],q[2],q[3]},
  Rmatrix3D{
  {1-2*(q[2]*q[2]+q[3]*q[3]), 2*(q[1]*q[2]-q[3]*q[0]), 2*(q[1]*q[3]+q[2]*q[0])},
  {2*(q[1]*q[2]+q[3]*q[0]), 1-2*(q[1]*q[1]+q[3]*q[3]), 2*(q[2]*q[3]-q[1]*q[0])},
  {2*(q[1]*q[3]-q[2]*q[0]), 2*(q[2]*q[3]+q[1]*q[0]), 1-2*(q[1]*q[1]+q[2]*q[2])}
  } { }

  virtual ~PutFishOnBlocks() {}

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

  void operator()(const Real , const Real , const Real , const Real,
                  ObstacleBlock*const,
                  const std::vector<VolumeSegment_OBB*>&) const;
  virtual void constructSurface(const Real , const Real , const Real , const Real,
                  ObstacleBlock*const,
                  const std::vector<VolumeSegment_OBB*>&) const;
  virtual void constructInternl(const Real , const Real , const Real , const Real,
                  ObstacleBlock*const,
                  const std::vector<VolumeSegment_OBB*>&) const;
  virtual void signedDistanceSqrt(ObstacleBlock*const) const;
};

struct PutNacaOnBlocks: public PutFishOnBlocks
{
  PutNacaOnBlocks(FishMidlineData* _cfish, const Real p[3], const Real q[4]): PutFishOnBlocks(_cfish, p, q) { }

  Real getSmallerDistToMidLPlanar(const int start_s, const Real x[3], int & final_s) const;

  void constructSurface(const Real , const Real , const Real , const Real,
                  ObstacleBlock*const,
                  const std::vector<VolumeSegment_OBB*>&) const override;
  void constructInternl(const Real , const Real , const Real , const Real,
                  ObstacleBlock*const,
                  const std::vector<VolumeSegment_OBB*>&) const override;
};

CubismUP_3D_NAMESPACE_END
