//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include "FishLibrary.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];
using CHIMAT = Real [CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];


void FishMidlineData::integrateLinearMomentum()
{
  // Compute the center of mass and center of mass velocities of the fish.
  // Then change midline frame of reference to that (new origin is the fish center of mass).
  // Mass is computed as:
  // M = int_{x} int_{y} int_{z} dx dy dz = int_{s} int_{E} |Jacobian| ds dE
  // where E is an elliptic cross-section of the fish.
  // The coordinate transformation that gives us the Jacobian is:
  // x(s,h1,h2) = rS(s) + nor(s)*width(s)*h1 + bin(s)*height(s)*h2
  // where h1,h2 are the coordinates in the ellipse
  // Center of mass (and its velocity) is computed as:
  // C_{x} = 1/M * int_{x} int_{y} int_{z} x dx dy dz = int_{s} int_{E} x |Jacobian| ds dE

  Real V=0, cmx=0, cmy=0, cmz=0, lmx=0, lmy=0, lmz=0;
  #pragma omp parallel for schedule(static) reduction(+:V,cmx,cmy,cmz,lmx,lmy,lmz)
  for(int i=0;i<Nm;++i)
  {  
    const Real ds = 0.5* ( (i==0) ? rS[1]-rS[0] : ((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]) );

    const Real c0 = norY[i]*binZ[i]-norZ[i]*binY[i];
    const Real c1 = norZ[i]*binX[i]-norX[i]*binZ[i];
    const Real c2 = norX[i]*binY[i]-norY[i]*binX[i];

    const Real x0dot = _d_ds(i, rX  , Nm);
    const Real x1dot = _d_ds(i, rY  , Nm);
    const Real x2dot = _d_ds(i, rZ  , Nm);
    const Real n0dot = _d_ds(i, norX, Nm);
    const Real n1dot = _d_ds(i, norY, Nm);
    const Real n2dot = _d_ds(i, norZ, Nm);
    const Real b0dot = _d_ds(i, binX, Nm);
    const Real b1dot = _d_ds(i, binY, Nm);
    const Real b2dot = _d_ds(i, binZ, Nm);

    const Real w = width[i];
    const Real H = height[i];

    const Real aux1 = w*H         *(c0*x0dot+c1*x1dot+c2*x2dot)*ds;
    const Real aux2 = 0.25*w*w*w*H*(c0*n0dot+c1*n1dot+c2*n2dot)*ds;
    const Real aux3 = 0.25*w*H*H*H*(c0*b0dot+c1*b1dot+c2*b2dot)*ds;

    V   += aux1;
    cmx += rX[i]*aux1 +norX[i]*aux2+ binX[i]*aux3;
    cmy += rY[i]*aux1 +norY[i]*aux2+ binY[i]*aux3;
    cmz += rZ[i]*aux1 +norZ[i]*aux2+ binZ[i]*aux3;
    lmx += vX[i]*aux1+vNorX[i]*aux2+vBinX[i]*aux3;
    lmy += vY[i]*aux1+vNorY[i]*aux2+vBinY[i]*aux3;
    lmz += vZ[i]*aux1+vNorZ[i]*aux2+vBinZ[i]*aux3;
  }
  const Real volume = V*M_PI;
  const Real aux = M_PI/volume;
  cmx *= aux;
  cmy *= aux;
  cmz *= aux;
  lmx *= aux;
  lmy *= aux;
  lmz *= aux;

  #pragma omp parallel for schedule(static)
  for(int i=0;i<Nm;++i)
  {
    rX[i]-=cmx; rY[i]-=cmy; rZ[i]-=cmz;
    vX[i]-=lmx; vY[i]-=lmy; vZ[i]-=lmz;
  }
}

void FishMidlineData::integrateAngularMomentum(const Real dt)
{
  // Compute the moments of inertia and angular velocities of the fish.
  // See comments in FishMidlineData::integrateLinearMomentum.
  Real JXX = 0;
  Real JYY = 0;
  Real JZZ = 0;
  Real JXY = 0;
  Real JYZ = 0;
  Real JZX = 0;
  Real AM_X = 0;
  Real AM_Y = 0;
  Real AM_Z = 0;
  #pragma omp parallel for reduction (+:JXX,JYY,JZZ,JXY,JYZ,JZX,AM_X,AM_Y,AM_Z)
  for(int i=0;i<Nm;++i)
  {
    const Real ds = 0.5* ( (i==0) ? rS[1]-rS[0] : ((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]) );

    const Real c0 = norY[i]*binZ[i]-norZ[i]*binY[i];
    const Real c1 = norZ[i]*binX[i]-norX[i]*binZ[i];
    const Real c2 = norX[i]*binY[i]-norY[i]*binX[i];
    const Real x0dot = _d_ds(i, rX  , Nm);
    const Real x1dot = _d_ds(i, rY  , Nm);
    const Real x2dot = _d_ds(i, rZ  , Nm);
    const Real n0dot = _d_ds(i, norX, Nm);
    const Real n1dot = _d_ds(i, norY, Nm);
    const Real n2dot = _d_ds(i, norZ, Nm);
    const Real b0dot = _d_ds(i, binX, Nm);
    const Real b1dot = _d_ds(i, binY, Nm);
    const Real b2dot = _d_ds(i, binZ, Nm);

    const Real M00 = width[i]*height[i];
    const Real M11 = 0.25*width[i]* width[i]* width[i]*height[i];
    const Real M22 = 0.25*width[i]*height[i]*height[i]*height[i];

    const Real cR = c0*x0dot + c1*x1dot + c2*x2dot;
    const Real cN = c0*n0dot + c1*n1dot + c2*n2dot;
    const Real cB = c0*b0dot + c1*b1dot + c2*b2dot;

    JXY += -ds*(cR* (rX[i]*rY[i]*M00 + norX[i]*norY[i]*M11 + binX[i]*binY[i]*M22) + cN*M11*(rX[i]*norY[i] + rY[i]*norX[i]) + cB*M22*(rX[i]*binY[i] + rY[i]*binX[i]));
    JZX += -ds*(cR* (rZ[i]*rX[i]*M00 + norZ[i]*norX[i]*M11 + binZ[i]*binX[i]*M22) + cN*M11*(rZ[i]*norX[i] + rX[i]*norZ[i]) + cB*M22*(rZ[i]*binX[i] + rX[i]*binZ[i]));
    JYZ += -ds*(cR* (rY[i]*rZ[i]*M00 + norY[i]*norZ[i]*M11 + binY[i]*binZ[i]*M22) + cN*M11*(rY[i]*norZ[i] + rZ[i]*norY[i]) + cB*M22*(rY[i]*binZ[i] + rZ[i]*binY[i]));

    const Real XX = ds*(cR* (rX[i]*rX[i]*M00 + norX[i]*norX[i]*M11 + binX[i]*binX[i]*M22) + cN*M11* (rX[i]*norX[i] + rX[i]*norX[i]) + cB*M22* (rX[i]*binX[i] + rX[i]*binX[i]));
    const Real YY = ds*(cR* (rY[i]*rY[i]*M00 + norY[i]*norY[i]*M11 + binY[i]*binY[i]*M22) + cN*M11* (rY[i]*norY[i] + rY[i]*norY[i]) + cB*M22* (rY[i]*binY[i] + rY[i]*binY[i]));
    const Real ZZ = ds*(cR* (rZ[i]*rZ[i]*M00 + norZ[i]*norZ[i]*M11 + binZ[i]*binZ[i]*M22) + cN*M11* (rZ[i]*norZ[i] + rZ[i]*norZ[i]) + cB*M22* (rZ[i]*binZ[i] + rZ[i]*binZ[i]));
    JXX += YY + ZZ;
    JYY += ZZ + XX;
    JZZ += YY + XX;

    const Real xd_y = cR* (vX[i]*rY[i]*M00 + vNorX[i]*norY[i]*M11 + vBinX[i]*binY[i]*M22) + cN*M11* (vX[i]*norY[i] + rY[i]*vNorX[i]) + cB*M22* (vX[i]*binY[i] + rY[i]*vBinX[i]);
    const Real x_yd = cR* (rX[i]*vY[i]*M00 + norX[i]*vNorY[i]*M11 + binX[i]*vBinY[i]*M22) + cN*M11* (rX[i]*vNorY[i] + rY[i]*norX[i]) + cB*M22* (rX[i]*vBinY[i] + vY[i]*binX[i]);
    const Real xd_z = cR* (rZ[i]*vX[i]*M00 + norZ[i]*vNorX[i]*M11 + binZ[i]*vBinX[i]*M22) + cN*M11* (rZ[i]*vNorX[i] + vX[i]*norZ[i]) + cB*M22* (rZ[i]*vBinX[i] + vX[i]*binZ[i]);
    const Real x_zd = cR* (vZ[i]*rX[i]*M00 + vNorZ[i]*norX[i]*M11 + vBinZ[i]*binX[i]*M22) + cN*M11* (vZ[i]*norX[i] + rX[i]*vNorZ[i]) + cB*M22* (vZ[i]*binX[i] + rX[i]*vBinZ[i]);
    const Real yd_z = cR* (vY[i]*rZ[i]*M00 + vNorY[i]*norZ[i]*M11 + vBinY[i]*binZ[i]*M22) + cN*M11* (vY[i]*norZ[i] + rZ[i]*vNorY[i]) + cB*M22* (vY[i]*binZ[i] + rZ[i]*vBinY[i]);
    const Real y_zd = cR* (rY[i]*vZ[i]*M00 + norY[i]*vNorZ[i]*M11 + binY[i]*vBinZ[i]*M22) + cN*M11* (rY[i]*vNorZ[i] + vZ[i]*norY[i]) + cB*M22* (rY[i]*vBinZ[i] + vZ[i]*binY[i]);

    AM_X += (y_zd - yd_z)*ds;
    AM_Y += (xd_z - x_zd)*ds;
    AM_Z += (x_yd - xd_y)*ds;
  }

  const Real eps = std::numeric_limits<Real>::epsilon();
  if (JXX < eps) JXX +=eps;
  if (JYY < eps) JYY +=eps;
  if (JZZ < eps) JZZ +=eps;
  JXX  *= M_PI; JYY  *= M_PI; JZZ  *= M_PI;
  JXY  *= M_PI; JYZ  *= M_PI; JZX  *= M_PI;
  AM_X *= M_PI; AM_Y *= M_PI; AM_Z *= M_PI;

  //Invert I
  const Real m00 = JXX; const Real m01 = JXY; const Real m02 = JZX;
  const Real m11 = JYY; const Real m12 = JYZ; const Real m22 = JZZ;
  const Real a00 = m22*m11 - m12*m12;
  const Real a01 = m02*m12 - m22*m01;
  const Real a02 = m01*m12 - m02*m11;
  const Real a11 = m22*m00 - m02*m02;
  const Real a12 = m01*m02 - m00*m12;
  const Real a22 = m00*m11 - m01*m01;
  const Real determinant =  1.0/((m00 * a00) + (m01 * a01) + (m02 * a02));
  
  angvel_internal[0] = (a00*AM_X + a01*AM_Y + a02*AM_Z)*determinant;
  angvel_internal[1] = (a01*AM_X + a11*AM_Y + a12*AM_Z)*determinant;
  angvel_internal[2] = (a02*AM_X + a12*AM_Y + a22*AM_Z)*determinant;
  const Real dqdt[4] = {
    0.5*( - angvel_internal[0]*quaternion_internal[1] - angvel_internal[1]*quaternion_internal[2] - angvel_internal[2]*quaternion_internal[3] ),
    0.5*( + angvel_internal[0]*quaternion_internal[0] + angvel_internal[1]*quaternion_internal[3] - angvel_internal[2]*quaternion_internal[2] ),
    0.5*( - angvel_internal[0]*quaternion_internal[3] + angvel_internal[1]*quaternion_internal[0] + angvel_internal[2]*quaternion_internal[1] ),
    0.5*( + angvel_internal[0]*quaternion_internal[2] - angvel_internal[1]*quaternion_internal[1] + angvel_internal[2]*quaternion_internal[0] )
  };
  quaternion_internal[0] -= dt * dqdt[0];
  quaternion_internal[1] -= dt * dqdt[1];
  quaternion_internal[2] -= dt * dqdt[2];
  quaternion_internal[3] -= dt * dqdt[3];
  const Real invD = 1.0/std::sqrt(quaternion_internal[0]*quaternion_internal[0] + quaternion_internal[1]*quaternion_internal[1] 
                                + quaternion_internal[2]*quaternion_internal[2] + quaternion_internal[3]*quaternion_internal[3]);
  quaternion_internal[0] *= invD;
  quaternion_internal[1] *= invD;
  quaternion_internal[2] *= invD;
  quaternion_internal[3] *= invD;

  //now we do the rotation
  Real R[3][3];
  R[0][0] = 1-2*(quaternion_internal[2]*quaternion_internal[2]+quaternion_internal[3]*quaternion_internal[3]);
  R[0][1] =   2*(quaternion_internal[1]*quaternion_internal[2]-quaternion_internal[3]*quaternion_internal[0]);
  R[0][2] =   2*(quaternion_internal[1]*quaternion_internal[3]+quaternion_internal[2]*quaternion_internal[0]);

  R[1][0] =   2*(quaternion_internal[1]*quaternion_internal[2]+quaternion_internal[3]*quaternion_internal[0]); 
  R[1][1] = 1-2*(quaternion_internal[1]*quaternion_internal[1]+quaternion_internal[3]*quaternion_internal[3]);
  R[1][2] =   2*(quaternion_internal[2]*quaternion_internal[3]-quaternion_internal[1]*quaternion_internal[0]);

  R[2][0] =   2*(quaternion_internal[1]*quaternion_internal[3]-quaternion_internal[2]*quaternion_internal[0]); 
  R[2][1] =   2*(quaternion_internal[2]*quaternion_internal[3]+quaternion_internal[1]*quaternion_internal[0]); 
  R[2][2] = 1-2*(quaternion_internal[1]*quaternion_internal[1]+quaternion_internal[2]*quaternion_internal[2]);
  #pragma omp parallel for schedule(static)
  for(int i=0;i<Nm;++i)
  {
    //rotation position and velocity
    {
      Real p[3] = {rX[i],rY[i],rZ[i]};
      rX[i] = R[0][0] * p[0]  + R[0][1] * p[1] + R[0][2] * p[2];
      rY[i] = R[1][0] * p[0]  + R[1][1] * p[1] + R[1][2] * p[2];
      rZ[i] = R[2][0] * p[0]  + R[2][1] * p[1] + R[2][2] * p[2];
      Real v[3] = {vX[i],vY[i],vZ[i]};
      vX[i] = R[0][0] * v[0]  + R[0][1] * v[1] + R[0][2] * v[2];
      vY[i] = R[1][0] * v[0]  + R[1][1] * v[1] + R[1][2] * v[2];
      vZ[i] = R[2][0] * v[0]  + R[2][1] * v[1] + R[2][2] * v[2];
      vX[i] +=  angvel_internal[2]*rY[i] - angvel_internal[1]*rZ[i];
      vY[i] +=  angvel_internal[0]*rZ[i] - angvel_internal[2]*rX[i];
      vZ[i] +=  angvel_internal[1]*rX[i] - angvel_internal[0]*rY[i];
    }
    //rotation normal vector
    {
      Real p[3] = {norX[i],norY[i],norZ[i]};
      norX[i] = R[0][0] * p[0]  + R[0][1] * p[1] + R[0][2] * p[2];
      norY[i] = R[1][0] * p[0]  + R[1][1] * p[1] + R[1][2] * p[2];
      norZ[i] = R[2][0] * p[0]  + R[2][1] * p[1] + R[2][2] * p[2];
      Real v[3] = {vNorX[i],vNorY[i],vNorZ[i]};
      vNorX[i] = R[0][0] * v[0]  + R[0][1] * v[1] + R[0][2] * v[2];
      vNorY[i] = R[1][0] * v[0]  + R[1][1] * v[1] + R[1][2] * v[2];
      vNorZ[i] = R[2][0] * v[0]  + R[2][1] * v[1] + R[2][2] * v[2];
      vNorX[i] += angvel_internal[2]*norY[i] - angvel_internal[1]*norZ[i];
      vNorY[i] += angvel_internal[0]*norZ[i] - angvel_internal[2]*norX[i];
      vNorZ[i] += angvel_internal[1]*norX[i] - angvel_internal[0]*norY[i];
    }
    //rotation binormal vector
    {
      Real p[3] = {binX[i],binY[i],binZ[i]};
      binX[i] = R[0][0] * p[0]  + R[0][1] * p[1] + R[0][2] * p[2];
      binY[i] = R[1][0] * p[0]  + R[1][1] * p[1] + R[1][2] * p[2];
      binZ[i] = R[2][0] * p[0]  + R[2][1] * p[1] + R[2][2] * p[2];
      Real v[3] = {vBinX[i],vBinY[i],vBinZ[i]};
      vBinX[i] = R[0][0] * v[0]  + R[0][1] * v[1] + R[0][2] * v[2];
      vBinY[i] = R[1][0] * v[0]  + R[1][1] * v[1] + R[1][2] * v[2];
      vBinZ[i] = R[2][0] * v[0]  + R[2][1] * v[1] + R[2][2] * v[2];
      vBinX[i] += angvel_internal[2]*binY[i] - angvel_internal[1]*binZ[i];
      vBinY[i] += angvel_internal[0]*binZ[i] - angvel_internal[2]*binX[i];
      vBinZ[i] += angvel_internal[1]*binX[i] - angvel_internal[0]*binY[i];
    }
  }
}

void VolumeSegment_OBB::prepare(std::pair<int, int> _s_range, const Real bbox[3][2], const Real h)
{
  safe_distance = (SURFDH+2)*h; //two points on each side for Towers
  s_range.first = _s_range.first;
  s_range.second = _s_range.second;
  for(int i=0; i<3; ++i) {
    w[i] = (bbox[i][1]-bbox[i][0])/2 + safe_distance;
    c[i] = (bbox[i][1]+bbox[i][0])/2;
    assert(w[i]>0);
  }
}

void VolumeSegment_OBB::normalizeNormals()
{
  const Real magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
  const Real magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
  const Real magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
  assert(magI > std::numeric_limits<Real>::epsilon());
  assert(magJ > std::numeric_limits<Real>::epsilon());
  assert(magK > std::numeric_limits<Real>::epsilon());
  const Real invMagI = Real(1)/magI;
  const Real invMagJ = Real(1)/magJ;
  const Real invMagK = Real(1)/magK;

  for(int i=0;i<3;++i) {
    // also take absolute value since thats what we need when doing intersection checks later
    normalI[i]=std::fabs(normalI[i])*invMagI;
    normalJ[i]=std::fabs(normalJ[i])*invMagJ;
    normalK[i]=std::fabs(normalK[i])*invMagK;
  }
}

void VolumeSegment_OBB::changeToComputationalFrame(const Real position[3], const Real quaternion[4])
{
  // we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
  const Real a = quaternion[0];
  const Real x = quaternion[1];
  const Real y = quaternion[2];
  const Real z = quaternion[3];
  const Real Rmatrix[3][3] = {
      {(Real)1.-2*(y*y+z*z),(Real)    2*(x*y-z*a),(Real)    2*(x*z+y*a)},
      {(Real)   2*(x*y+z*a),(Real) 1.-2*(x*x+z*z),(Real)    2*(y*z-x*a)},
      {(Real)   2*(x*z-y*a),(Real)    2*(y*z+x*a),(Real) 1.-2*(x*x+y*y)}
  };
  const Real p[3] = {c[0],c[1],c[2]};
  const Real nx[3] = {normalI[0],normalI[1],normalI[2]};
  const Real ny[3] = {normalJ[0],normalJ[1],normalJ[2]};
  const Real nz[3] = {normalK[0],normalK[1],normalK[2]};
  for(int i=0;i<3;++i) {
    c[i]      = Rmatrix[i][0]*p[0]  +Rmatrix[i][1]*p[1]  +Rmatrix[i][2]*p[2];
    normalI[i]= Rmatrix[i][0]*nx[0] +Rmatrix[i][1]*nx[1] +Rmatrix[i][2]*nx[2];
    normalJ[i]= Rmatrix[i][0]*ny[0] +Rmatrix[i][1]*ny[1] +Rmatrix[i][2]*ny[2];
    normalK[i]= Rmatrix[i][0]*nz[0] +Rmatrix[i][1]*nz[1] +Rmatrix[i][2]*nz[2];
  }
  c[0] +=position[0];
  c[1] +=position[1];
  c[2] +=position[2];

  normalizeNormals();
  assert(normalI[0]>=0 && normalI[1]>=0 && normalI[2]>=0);
  assert(normalJ[0]>=0 && normalJ[1]>=0 && normalJ[2]>=0);
  assert(normalK[0]>=0 && normalK[1]>=0 && normalK[2]>=0);

  // Find the x,y,z max extents in lab frame ( exploit normal(I,J,K)[:] >=0 )
  const Real widthXvec[] = {w[0]*normalI[0], w[0]*normalI[1], w[0]*normalI[2]};
  const Real widthYvec[] = {w[1]*normalJ[0], w[1]*normalJ[1], w[1]*normalJ[2]};
  const Real widthZvec[] = {w[2]*normalK[0], w[2]*normalK[1], w[2]*normalK[2]};

  for(int i=0; i<3; ++i) {
    objBoxLabFr[i][0] = c[i] -widthXvec[i] -widthYvec[i] -widthZvec[i];
    objBoxLabFr[i][1] = c[i] +widthXvec[i] +widthYvec[i] +widthZvec[i];
    objBoxObjFr[i][0] = c[i] -w[i];
    objBoxObjFr[i][1] = c[i] +w[i];
  }
}

#define DBLCHECK
bool VolumeSegment_OBB::isIntersectingWithAABB(const Real start[3],const Real end[3]) const
{
  // Remember: Incoming coordinates are cell centers, not cell faces
  //start and end are two diagonally opposed corners of grid block
  // GN halved the safety here but added it back to w[] in prepare
  const Real AABB_w[3] = { //half block width + safe distance
      (end[0] - start[0])/2 + safe_distance,
      (end[1] - start[1])/2 + safe_distance,
      (end[2] - start[2])/2 + safe_distance
  };

  const Real AABB_c[3] = { //block center
    (end[0] + start[0])/2,
    (end[1] + start[1])/2,
    (end[2] + start[2])/2
  };

  const Real AABB_box[3][2] = {
    {AABB_c[0] - AABB_w[0],  AABB_c[0] + AABB_w[0]},
    {AABB_c[1] - AABB_w[1],  AABB_c[1] + AABB_w[1]},
    {AABB_c[2] - AABB_w[2],  AABB_c[2] + AABB_w[2]}
  };

  assert(AABB_w[0]>0 && AABB_w[1]>0 && AABB_w[2]>0);

  // Now Identify the ones that do not intersect
  using std::max; using std::min;
  Real intersectionLabFrame[3][2] = {
  {max(objBoxLabFr[0][0],AABB_box[0][0]),min(objBoxLabFr[0][1],AABB_box[0][1])},
  {max(objBoxLabFr[1][0],AABB_box[1][0]),min(objBoxLabFr[1][1],AABB_box[1][1])},
  {max(objBoxLabFr[2][0],AABB_box[2][0]),min(objBoxLabFr[2][1],AABB_box[2][1])}
  };

  if ( intersectionLabFrame[0][1] - intersectionLabFrame[0][0] < 0
    || intersectionLabFrame[1][1] - intersectionLabFrame[1][0] < 0
    || intersectionLabFrame[2][1] - intersectionLabFrame[2][0] < 0 )
    return false;

  #ifdef DBLCHECK
    const Real widthXbox[3] = {AABB_w[0]*normalI[0], AABB_w[0]*normalJ[0], AABB_w[0]*normalK[0]}; // This is x-width of box, expressed in fish frame
    const Real widthYbox[3] = {AABB_w[1]*normalI[1], AABB_w[1]*normalJ[1], AABB_w[1]*normalK[1]}; // This is y-width of box, expressed in fish frame
    const Real widthZbox[3] = {AABB_w[2]*normalI[2], AABB_w[2]*normalJ[2], AABB_w[2]*normalK[2]}; // This is z-height of box, expressed in fish frame

    const Real boxBox[3][2] = {
      { AABB_c[0] -widthXbox[0] -widthYbox[0] -widthZbox[0],
        AABB_c[0] +widthXbox[0] +widthYbox[0] +widthZbox[0]},
      { AABB_c[1] -widthXbox[1] -widthYbox[1] -widthZbox[1],
        AABB_c[1] +widthXbox[1] +widthYbox[1] +widthZbox[1]},
      { AABB_c[2] -widthXbox[2] -widthYbox[2] -widthZbox[2],
        AABB_c[2] +widthXbox[2] +widthYbox[2] +widthZbox[2]}
    };

    Real intersectionFishFrame[3][2] = {
     {max(boxBox[0][0],objBoxObjFr[0][0]), min(boxBox[0][1],objBoxObjFr[0][1])},
     {max(boxBox[1][0],objBoxObjFr[1][0]), min(boxBox[1][1],objBoxObjFr[1][1])},
     {max(boxBox[2][0],objBoxObjFr[2][0]), min(boxBox[2][1],objBoxObjFr[2][1])}
    };

    if ( intersectionFishFrame[0][1] - intersectionFishFrame[0][0] < 0
      || intersectionFishFrame[1][1] - intersectionFishFrame[1][0] < 0
      || intersectionFishFrame[2][1] - intersectionFishFrame[2][0] < 0 )
      return false;
  #endif

  return true;
}

void PutFishOnBlocks::operator()(const Real h, const Real ox, const Real oy, const Real oz, ObstacleBlock* const oblock, const std::vector<VolumeSegment_OBB*>& vSegments) const
{
  auto & sdf_array = oblock->sdfLab;
  for(int iz=0; iz<ScalarBlock::sizeZ+2; ++iz)
  for(int iy=0; iy<ScalarBlock::sizeY+2; ++iy)
  for(int ix=0; ix<ScalarBlock::sizeX+2; ++ix)
    sdf_array[iz][iy][ix] = -1;
  constructSurface  (h, ox, oy, oz, oblock, vSegments);
  constructInternl  (h, ox, oy, oz, oblock, vSegments);
  signedDistanceSqrt(oblock);
}

inline Real distPlane(const Real p1[3], const Real p2[3], const Real p3[3],
                      const Real s[3], const Real IN[3])
{
  // make p1 origin of a frame of ref
  const Real t[3] = {  s[0]-p1[0],  s[1]-p1[1],  s[2]-p1[2] };
  const Real u[3] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
  const Real v[3] = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };
  const Real i[3] = { IN[0]-p1[0], IN[1]-p1[1], IN[2]-p1[2] };
  // normal to the plane:
  const Real n[3] = {  u[1]*v[2] - u[2]*v[1],
                       u[2]*v[0] - u[0]*v[2],
                       u[0]*v[1] - u[1]*v[0]};
  // if normal points inside then this is going to be positive:
  const Real projInner = i[0]*n[0] + i[1]*n[1] + i[2]*n[2];
  // if normal points outside we need to change sign of result:
  const Real signIn = projInner>0 ? 1 : -1;
  //every point of the plane will have no projection onto n
  // therefore, distance of t from plane is:
  const Real norm = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  return signIn * (t[0]*n[0] + t[1]*n[1] + t[2]*n[2]) / norm;
}

void PutFishOnBlocks::constructSurface(const Real h, const Real ox, const Real oy, const Real oz, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB*>&vSegments) const
{
  //Construct the surface of the fish.
  //By definition, this is where SDF = 0.
  //Here, we know an analytical expression for the fish surface:
  // x(s,theta) = r(s) + width(s)*cos(theta)*nor(s) + height(s)*sin(theta)*bin(s)
  //We loop over the discretized version of this equation and for each point of the surface we find
  //the grid point that is the closest.
  //The SDF is computed for that grid point plus the grid points that are within
  //a distance of 2h of that point. This is done because we need the SDF in a zone of +-2h of the 
  //actual surface, so that we can use the mollified Heaviside function from Towers for chi.

  //Pointers to discrete values that describe the surface.
  const Real *const rX = cfish->rX, *const norX  = cfish->norX , *const vBinX = cfish->vBinX;
  const Real *const rY = cfish->rY, *const norY  = cfish->norY , *const vBinY = cfish->vBinY;
  const Real *const rZ = cfish->rZ, *const norZ  = cfish->norZ , *const vBinZ = cfish->vBinZ;
  const Real *const vX = cfish->vX, *const vNorX = cfish->vNorX, *const binX  = cfish->binX ;
  const Real *const vY = cfish->vY, *const vNorY = cfish->vNorY, *const binY  = cfish->binY ;
  const Real *const vZ = cfish->vZ, *const vNorZ = cfish->vNorZ, *const binZ  = cfish->binZ ;
  Real *const width = cfish->width;
  Real *const height = cfish->height;

  //These are typically 8x8x8 (blocksize=8) matrices that are filled here.
  CHIMAT & __restrict__ CHI = defblock->chi;
  UDEFMAT & __restrict__ UDEF = defblock->udef;

  //This is an (8+2)x(8+2)x(8+2) matrix with the SDF. We compute the sdf for the 8x8x8 block
  //but we also add +-1 grid points of ghost values. That way there is no need for communication
  //later on. We do this because an analytical expression is available and the extra cost
  //for this computation is small.
  auto & __restrict__ SDFLAB = defblock->sdfLab;

  //Origin of the block, displaced by (-h,-h,-h). This is where the first ghost cell is.
  const Real org[3] = {ox-h,oy-h,oz-h};

  const Real invh = 1.0/h;
  const int BS[3] = {ScalarBlock::sizeX+2, ScalarBlock::sizeY+2, ScalarBlock::sizeZ+2};

  // needed variables to store sensor location
  const Real *const rS = cfish->rS;
  const Real length = cfish->length;

  // save location for tip of head
  Real myP[3] ={ rX[0], rY[0], rZ[0] };
  changeToComputationalFrame(myP);
  cfish->sensorLocation[0] = myP[0];
  cfish->sensorLocation[1] = myP[1];
  cfish->sensorLocation[2] = myP[2];

  //Loop over vSegments of this block. 
  for(size_t i=0; i<vSegments.size(); ++i)
  {
    //Each segment has a range of s for the surface that intersects this block.
    const int firstSegm = std::max(vSegments[i]->s_range.first,            1);
    const int lastSegm =  std::min(vSegments[i]->s_range.second, cfish->Nm-2);

    //Loop over discrete s of surface x(s,theta)
    for(int ss=firstSegm; ss<=lastSegm; ++ss) 
    {
      if (height[ss]<=0) height[ss]=1e-10;
      if (width [ss]<=0) width [ss]=1e-10;

      //Here we discretize theta for the given s (for given s, the cross-section is an ellipse) 
      //This is done by setting the maximum arclength to h/2.
      //Note that maximum arclength = sin(dtheta) * (major_axis+h) = h/2
      const Real major_axis = std::max(height[ss], width[ss]);
      const Real dtheta_tgt = std::fabs(std::asin(h/(major_axis+h)/2));
      int Ntheta = std::ceil(2*M_PI/dtheta_tgt);
      if (Ntheta % 2 == 1) Ntheta++; //force Ntheta to be even so that the fish is symmetric
      const Real dtheta = 2*M_PI/((Real) Ntheta);

      //theta = 0 at the major axis, this variable takes care of that
      const Real offset = height[ss] > width[ss] ? M_PI/2 : 0;

      //Loop over discrete theta of surface x(s,theta)
      for(int tt=0; tt<Ntheta; ++tt)
      {
        const Real theta = tt*dtheta + offset;
        const Real sinth = std::sin(theta), costh = std::cos(theta);

        //Current surface point (in frame of reference of fish)
        myP[0] = rX[ss] +width[ss]*costh*norX[ss]+height[ss]*sinth*binX[ss];
        myP[1] = rY[ss] +width[ss]*costh*norY[ss]+height[ss]*sinth*binY[ss];
        myP[2] = rZ[ss] +width[ss]*costh*norZ[ss]+height[ss]*sinth*binZ[ss];

        changeToComputationalFrame(myP);

        // save location for side of head; for angle = 0 and angle = pi this is a sensor location
        if( rS[ss] <= 0.04*length && rS[ss+1] > 0.04*length )
        {
          if( tt == 0 )
          {
            cfish->sensorLocation[1*3+0] = myP[0];
            cfish->sensorLocation[1*3+1] = myP[1];
            cfish->sensorLocation[1*3+2] = myP[2];
          }
          if( tt == (int)Ntheta/2 )
          {
            cfish->sensorLocation[2*3+0] = myP[0];
            cfish->sensorLocation[2*3+1] = myP[1];
            cfish->sensorLocation[2*3+2] = myP[2];
          }
        }

        // Find index of nearest grid point to myP
        const int iap[3] = {(int)std::floor((myP[0]-org[0])*invh),
                            (int)std::floor((myP[1]-org[1])*invh),
                            (int)std::floor((myP[2]-org[2])*invh)};

        // Loop over that point and a neighborhood of +-3 points, to compute the SDF near the surface
        const int nei = 3;
        const int ST[3] = { iap[0]-nei , iap[1]-nei, iap[2]-nei};
        const int EN[3] = { iap[0]+nei , iap[1]+nei, iap[2]+nei};

        if(EN[0] <= 0 || ST[0] > BS[0]) continue; // NearNeigh loop
        if(EN[1] <= 0 || ST[1] > BS[1]) continue; // does not intersect
        if(EN[2] <= 0 || ST[2] > BS[2]) continue; // with this block

        //Store the surface point at the next and the previous cross-sections.
        //They will be used to compute the SDF later.
        Real pP[3] = {rX[ss+1] +width[ss+1]*costh*norX[ss+1]+height[ss+1]*sinth*binX[ss+1],
                      rY[ss+1] +width[ss+1]*costh*norY[ss+1]+height[ss+1]*sinth*binY[ss+1],
                      rZ[ss+1] +width[ss+1]*costh*norZ[ss+1]+height[ss+1]*sinth*binZ[ss+1]};
        Real pM[3] = {rX[ss-1] +width[ss-1]*costh*norX[ss-1]+height[ss-1]*sinth*binX[ss-1],
                      rY[ss-1] +width[ss-1]*costh*norY[ss-1]+height[ss-1]*sinth*binY[ss-1],
                      rZ[ss-1] +width[ss-1]*costh*norZ[ss-1]+height[ss-1]*sinth*binZ[ss-1]};
        changeToComputationalFrame(pM);
        changeToComputationalFrame(pP);

        //Deformation velocity of surface point.
        Real udef[3] = { vX[ss] +width[ss]*costh*vNorX[ss]+height[ss]*sinth*vBinX[ss],
                         vY[ss] +width[ss]*costh*vNorY[ss]+height[ss]*sinth*vBinY[ss],
                         vZ[ss] +width[ss]*costh*vNorZ[ss]+height[ss]*sinth*vBinZ[ss]};
        changeVelocityToComputationalFrame(udef);

        for(int sz = std::max(0, ST[2]); sz < std::min(EN[2], BS[2]); ++sz)
        for(int sy = std::max(0, ST[1]); sy < std::min(EN[1], BS[1]); ++sy)
        for(int sx = std::max(0, ST[0]); sx < std::min(EN[0], BS[0]); ++sx)
        {
          //Grid point in the neighborhood near surface
          Real p[3];
          p[0] = ox + h * (sx - 1 + 0.5);
          p[1] = oy + h * (sy - 1 + 0.5);
          p[2] = oz + h * (sz - 1 + 0.5);

          const Real dist0 = eulerDistSq3D(p, myP);
          const Real distP = eulerDistSq3D(p, pP);
          const Real distM = eulerDistSq3D(p, pM);

          // check if this grid point has already found a closer surf-point:
          if(std::fabs(SDFLAB[sz][sy][sx])<std::min({dist0,distP,distM})) continue;

          // if this grid point is > 2h distance of the grid point that is nearest to the surface
          // don't compute the sdf
          if(std::min({dist0,distP,distM})>4*h*h) continue;

          changeFromComputationalFrame(p);

          // among the three points myP, pP, pM find the two that are the closest to this grid point
          int close_s = ss, secnd_s = ss + (distP<distM? 1 : -1);
          Real dist1 = dist0, dist2 = distP<distM? distP : distM;
          if(distP < dist0 || distM < dist0)
          {
            dist1 = dist2; 
            dist2 = dist0;
            close_s = secnd_s;
            secnd_s = ss;
          }

          // Interpolate the surface deformation velocity to this grid point.
          // W behaves like hat interpolation kernel that is used for internal
          // fish points. Introducing W (used to be W=1) smoothens transition
          // from surface to internal points. In fact, later we plus equal
          // udef*hat of internal points. If hat>0, point should behave like
          // internal point, meaning that fish-section udef rotation should
          // multiply distance from midline instead of entire half-width.
          // Remember that uder will become udef / chi, so W simplifies out.
          const Real W = std::max(1 - std::sqrt(dist1) * (invh / 3), (Real)0);
          const bool inRange =(sz - 1 >=0 && sz - 1 < ScalarBlock::sizeZ &&
                               sy - 1 >=0 && sy - 1 < ScalarBlock::sizeY &&
                               sx - 1 >=0 && sx - 1 < ScalarBlock::sizeX);
          if (inRange)
          {
            UDEF[sz-1][sy-1][sx-1][0] = W * udef[0];
            UDEF[sz-1][sy-1][sx-1][1] = W * udef[1];
            UDEF[sz-1][sy-1][sx-1][2] = W * udef[2];
            CHI [sz-1][sy-1][sx-1] = W; // Not chi, just used to interpolate udef!
          }

          //Now we compute the SDF of that point.
          //If that point is close to the tail, we project is onto the plane defined by the tail
          //and then compute the sign of the distance based on the side of the plane the point lies on.

          //Else, we model the span between two ellipses (cross-sections) as a spherical segment.
          //See also: http://mathworld.wolfram.com/SphericalSegment.html
          //The spherical segment is defined by two circles.
          //To define those cyrcles, we need their centers and a point on them.
          //The points myP,pM,pP (we are using two of those three here) are the cyrcle points we need.
          //The centers are found by taking the normal starting from myP/pM/pP to the line defined
          //by the vector R1:
          const Real R1[3]= {rX[secnd_s]-rX[close_s], rY[secnd_s]-rY[close_s], rZ[secnd_s]-rZ[close_s]};
          const Real normR1 = 1.0/(1e-21+std::sqrt(R1[0]*R1[0]+R1[1]*R1[1]+R1[2]*R1[2]));
          const Real nn[3] = {R1[0]*normR1,R1[1]*normR1,R1[2]*normR1};

          const Real P1[3]= {width[close_s]*costh*norX[close_s]+height[close_s]*sinth*binX[close_s],
                             width[close_s]*costh*norY[close_s]+height[close_s]*sinth*binY[close_s],
                             width[close_s]*costh*norZ[close_s]+height[close_s]*sinth*binZ[close_s]};
          const Real P2[3]= {width[secnd_s]*costh*norX[secnd_s]+height[secnd_s]*sinth*binX[secnd_s],
                             width[secnd_s]*costh*norY[secnd_s]+height[secnd_s]*sinth*binY[secnd_s],
                             width[secnd_s]*costh*norZ[secnd_s]+height[secnd_s]*sinth*binZ[secnd_s]};
          const Real dot1  = P1[0]*R1[0] + P1[1]*R1[1] + P1[2]*R1[2];
          const Real dot2  = P2[0]*R1[0] + P2[1]*R1[1] + P2[2]*R1[2];
          const Real base1 = dot1 * normR1;
          const Real base2 = dot2 * normR1;

          //These are a^2 and b^2 in http://mathworld.wolfram.com/SphericalSegment.html
          const Real radius_close  = std::pow( width[close_s]*costh,2) + std::pow(height[close_s]*sinth,2) - base1*base1;
          const Real radius_second = std::pow( width[secnd_s]*costh,2) + std::pow(height[secnd_s]*sinth,2) - base2*base2;

          const Real center_close  [3] = {rX[close_s]-nn[0]*base1,rY[close_s]-nn[1]*base1,rZ[close_s]-nn[2]*base1};
          const Real center_second [3] = {rX[secnd_s]+nn[0]*base2,rY[secnd_s]+nn[1]*base2,rZ[secnd_s]+nn[2]*base2};

          //This is h in http://mathworld.wolfram.com/SphericalSegment.html
          const Real dSsq = std::pow(center_close[0]-center_second[0], 2)
                           +std::pow(center_close[1]-center_second[1], 2)
                           +std::pow(center_close[2]-center_second[2], 2);

          const Real corr = 2*std::sqrt(radius_close*radius_second);

          if (close_s == cfish->Nm-2 || secnd_s == cfish->Nm-2) //point is close to tail
          {
            //compute the 5 corners of the pyramid around tail last point
            const int TT = cfish->Nm-1, TS = cfish->Nm-2;
            const Real PC[3] = {rX[TT], rY[TT], rZ[TT] };
            const Real PF[3] = {rX[TS], rY[TS], rZ[TS] };
            const Real DXT = p[0] - PF[0];
            const Real DYT = p[1] - PF[1];
            const Real DZT = p[2] - PF[2];
            const Real projW = ( width[TS]*norX[TS])*DXT+( width[TS]*norY[TS])*DYT+( width[TS]*norZ[TS])*DZT;
            const Real projH = (height[TS]*binX[TS])*DXT+(height[TS]*binY[TS])*DYT+(height[TS]*binZ[TS])*DZT;
            const int signW = projW > 0 ? 1 : -1;
            const int signH = projH > 0 ? 1 : -1;
            const Real PT[3] =  {rX[TS] + signH * height[TS]*binX[TS],
                                 rY[TS] + signH * height[TS]*binY[TS],
                                 rZ[TS] + signH * height[TS]*binZ[TS]};
            const Real PP[3] =  {rX[TS] + signW *  width[TS]*norX[TS],
                                 rY[TS] + signW *  width[TS]*norY[TS],
                                 rZ[TS] + signW *  width[TS]*norZ[TS]};
            SDFLAB[sz][sy][sx] = distPlane(PC, PT, PP, p, PF);
          }
          else if (dSsq >= radius_close+radius_second -corr) // if ds > delta radius
          { 
            // if the two cross-sections are close and have axis that do not differ much, we just
            // use the nearest neighbor to compute the sdf (no need for spherical segment model)
            const Real xMidl[3] = {rX[close_s], rY[close_s], rZ[close_s]};
            const Real grd2ML = eulerDistSq3D(p, xMidl);
            const Real sign = grd2ML > radius_close ? -1 : 1;
            SDFLAB[sz][sy][sx] = sign*dist1;
          }
          else //here we use the spherical segment model
          {
            const Real Rsq = (radius_close +radius_second -corr +dSsq) //radius of the spere
                            *(radius_close +radius_second +corr +dSsq)/4/dSsq;
            const Real maxAx = std::max(radius_close, radius_second);

            // 'submerged' fraction of radius:
            const Real d = std::sqrt((Rsq - maxAx)/dSsq); //(divided by ds)
            // position of the centre of the sphere:
            Real sign;
            if(radius_close> radius_second)
            {
              const Real xMidl[3] = {center_close[0] +(center_close[0]-center_second[0])*d,
                                     center_close[1] +(center_close[1]-center_second[1])*d,
                                     center_close[2] +(center_close[2]-center_second[2])*d};
              const Real grd2Core = eulerDistSq3D(p, xMidl);
              sign = grd2Core > Rsq ? -1 : 1;
            }
            else
            {
              const Real xMidl[3] = {center_second[0] +(center_second[0]-center_close[0])*d,
                                     center_second[1] +(center_second[1]-center_close[1])*d,
                                     center_second[2] +(center_second[2]-center_close[2])*d};
              const Real grd2Core = eulerDistSq3D(p, xMidl);
              sign = grd2Core > Rsq ? -1 : 1;
            }
            SDFLAB[sz][sy][sx] = sign*dist1;
          }
        }
      }
    }
  }
}

void PutFishOnBlocks::constructInternl(const Real h, const Real ox, const Real oy, const Real oz, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB*>&vSegments) const
{
  Real org[3] = {ox-h,oy-h,oz-h};
  const Real invh = 1.0/h;
  CHIMAT & __restrict__ CHI = defblock->chi;
  auto & __restrict__ SDFLAB = defblock->sdfLab;
  UDEFMAT & __restrict__ UDEF = defblock->udef;
  static constexpr int BS[3] = {ScalarBlock::sizeX+2, ScalarBlock::sizeY+2, ScalarBlock::sizeZ+2};
  const Real *const rX = cfish->rX, *const norX  = cfish->norX , *const vBinX = cfish->vBinX;
  const Real *const rY = cfish->rY, *const norY  = cfish->norY , *const vBinY = cfish->vBinY;
  const Real *const rZ = cfish->rZ, *const norZ  = cfish->norZ , *const vBinZ = cfish->vBinZ;
  const Real *const vX = cfish->vX, *const vNorX = cfish->vNorX, *const binX  = cfish->binX ;
  const Real *const vY = cfish->vY, *const vNorY = cfish->vNorY, *const binY  = cfish->binY ;
  const Real *const vZ = cfish->vZ, *const vNorZ = cfish->vNorZ, *const binZ  = cfish->binZ ;
  const Real *const width = cfish->width, *const height = cfish->height;

  // construct the deformation velocities (P2M with hat function as kernel)
  for(size_t i=0; i<vSegments.size(); ++i)
  {
  const int firstSegm = std::max(vSegments[i]->s_range.first,            1);
  const int lastSegm =  std::min(vSegments[i]->s_range.second, cfish->Nm-2);
  for(int ss=firstSegm; ss<=lastSegm; ++ss)
  {
    // P2M udef of a slice at this s
    const Real myWidth = width[ss], myHeight = height[ss];
    assert(myWidth > 0 && myHeight > 0);
    const int Nh = std::floor(myHeight/h); //floor bcz we already did interior
    for(int ih=-Nh+1; ih<Nh; ++ih)
    {
      const Real offsetH = ih * h;
      const Real currWidth = myWidth*std::sqrt(1-std::pow(offsetH/myHeight, 2));
      const int Nw = std::floor(currWidth/h);//floor bcz we already did interior
      for(int iw = -Nw+1; iw < Nw; ++iw)
      {
        const Real offsetW = iw * h;
        Real xp[3]= {rX[ss]+offsetW*norX[ss]+offsetH*binX[ss], 
                     rY[ss]+offsetW*norY[ss]+offsetH*binY[ss], 
                     rZ[ss]+offsetW*norZ[ss]+offsetH*binZ[ss]};
        changeToComputationalFrame(xp);
        xp[0] = (xp[0]-org[0])*invh; // how many grid points
        xp[1] = (xp[1]-org[1])*invh; // from this block origin
        xp[2] = (xp[2]-org[2])*invh; // is this fishpoint located at?
        const Real ap[3] = {
            std::floor(xp[0]), std::floor(xp[1]), std::floor(xp[2])
        };
        const int iap[3] = { (int)ap[0], (int)ap[1], (int)ap[2] };
        if(iap[0]+2 <= 0 || iap[0] >= BS[0]) continue; // hatP2M loop
        if(iap[1]+2 <= 0 || iap[1] >= BS[1]) continue; // does not intersect
        if(iap[2]+2 <= 0 || iap[2] >= BS[2]) continue; // with this block

        Real udef[3] = {vX[ss]+offsetW*vNorX[ss]+offsetH*vBinX[ss], 
                        vY[ss]+offsetW*vNorY[ss]+offsetH*vBinY[ss], 
                        vZ[ss]+offsetW*vNorZ[ss]+offsetH*vBinZ[ss]};
        changeVelocityToComputationalFrame(udef);
        Real wghts[3][2]; // P2M weights
        for(int c=0; c<3; ++c) {
          const Real t[2] = { // we floored, hat between xp and grid point +-1
              std::fabs(xp[c] -ap[c]), std::fabs(xp[c] -(ap[c] +1))
          };
          wghts[c][0] = 1.0 - t[0];
          wghts[c][1] = 1.0 - t[1];
          assert(wghts[c][0]>=0 && wghts[c][1]>=0);
        }

        for(int idz =std::max(0, iap[2]); idz <std::min(iap[2]+2, BS[2]); ++idz)
        for(int idy =std::max(0, iap[1]); idy <std::min(iap[1]+2, BS[1]); ++idy)
        for(int idx =std::max(0, iap[0]); idx <std::min(iap[0]+2, BS[0]); ++idx)
        {
          const int sx = idx - iap[0], sy = idy - iap[1], sz = idz - iap[2];
          assert( sx>=0 && sx<2 && sy>=0 && sy<2 && sz>=0 && sz<2 );
          const Real wxwywz = wghts[2][sz] * wghts[1][sy] * wghts[0][sx];
          assert(wxwywz>=0 && wxwywz<=1);

          if (idz -1 >=0 && idz-1 < ScalarBlock::sizeZ &&
              idy -1 >=0 && idy-1 < ScalarBlock::sizeY &&
              idx -1 >=0 && idx-1 < ScalarBlock::sizeX)
          {
               UDEF[idz-1][idy-1][idx-1][0] += wxwywz*udef[0];
               UDEF[idz-1][idy-1][idx-1][1] += wxwywz*udef[1];
               UDEF[idz-1][idy-1][idx-1][2] += wxwywz*udef[2];
               CHI [idz-1][idy-1][idx-1]    += wxwywz;
          }
          static constexpr Real eps = std::numeric_limits<Real>::epsilon();
          if( std::fabs(SDFLAB[idz][idy][idx]+1)<eps ) SDFLAB[idz][idy][idx] = 1;
          // set sign for all interior points
        }
      }
    }
  }
  }
}

void PutFishOnBlocks::signedDistanceSqrt(ObstacleBlock* const defblock) const
{
  static constexpr Real eps = std::numeric_limits<Real>::epsilon();
  auto & __restrict__ CHI    = defblock->chi;
  auto & __restrict__ UDEF   = defblock->udef;
  auto & __restrict__ SDFLAB = defblock->sdfLab;
  for(int iz=0; iz<ScalarBlock::sizeZ+2; iz++)
  for(int iy=0; iy<ScalarBlock::sizeY+2; iy++)
  for(int ix=0; ix<ScalarBlock::sizeX+2; ix++)
  {
    if (iz < ScalarBlock::sizeZ && iy < ScalarBlock::sizeY && ix < ScalarBlock::sizeX)
    {
      if (CHI[iz][iy][ix] > eps)
      {
        const Real normfac = 1.0/CHI[iz][iy][ix];
        UDEF[iz][iy][ix][0] *= normfac;
        UDEF[iz][iy][ix][1] *= normfac;
        UDEF[iz][iy][ix][2] *= normfac;
      } 
    }
    SDFLAB[iz][iy][ix] = SDFLAB[iz][iy][ix] >= 0 ? std::sqrt( SDFLAB[iz][iy][ix]) : -std::sqrt(-SDFLAB[iz][iy][ix]);
  }
}

void PutNacaOnBlocks::constructSurface(const Real h, const Real ox, const Real oy, const Real oz, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB*>&vSegments) const
{
  Real org[3] = {ox-h,oy-h,oz-h};
  const Real invh = 1.0/h;
  const Real* const rX = cfish->rX;
  const Real* const rY = cfish->rY;
  const Real* const norX = cfish->norX;
  const Real* const norY = cfish->norY;
  const Real* const vX = cfish->vX;
  const Real* const vY = cfish->vY;
  const Real* const vNorX = cfish->vNorX;
  const Real* const vNorY = cfish->vNorY;
  const Real* const width = cfish->width;
  const Real* const height = cfish->height;
  static constexpr int BS[3] = {ScalarBlock::sizeX+2, ScalarBlock::sizeY+2, ScalarBlock::sizeZ+2};
  CHIMAT & __restrict__ CHI = defblock->chi;
  auto & __restrict__ SDFLAB = defblock->sdfLab;
  UDEFMAT & __restrict__ UDEF = defblock->udef;

  // construct the shape (P2M with min(distance) as kernel) onto defblocks
  for(size_t i=0; i< vSegments.size(); ++i) {
    //iterate over segments contained in the vSegm intersecting this block:
    const int firstSegm = std::max(vSegments[i]->s_range.first,            1);
    const int lastSegm =  std::min(vSegments[i]->s_range.second, cfish->Nm-2);
    for(int ss=firstSegm; ss<=lastSegm; ++ss) {
      assert(height[ss]>0 && width[ss]>0);
      //for each segment, we have one point to left and right of midl
      for(int signp = -1; signp <= 1; signp+=2) {
        // create a surface point
        // special treatment of tail (width = 0 --> no ellipse, just line)
        Real myP[3] = {     rX[ss+0] +width[ss+0]*signp*norX[ss+0],
                            rY[ss+0] +width[ss+0]*signp*norY[ss+0], 0
        };
        const Real pP[3] = {rX[ss+1] +width[ss+1]*signp*norX[ss+1],
                            rY[ss+1] +width[ss+1]*signp*norY[ss+1], 0
        };
        const Real pM[3] = {rX[ss-1] +width[ss-1]*signp*norX[ss-1],
                            rY[ss-1] +width[ss-1]*signp*norY[ss-1], 0
        };
        changeToComputationalFrame(myP);
        const int iap[2] = {  (int)std::floor((myP[0]-org[0])*invh),
                              (int)std::floor((myP[1]-org[1])*invh)
        };
        Real udef[3] = { vX[ss+0] +width[ss+0]*signp*vNorX[ss+0],
                         vY[ss+0] +width[ss+0]*signp*vNorY[ss+0], 0
        };
        changeVelocityToComputationalFrame(udef);
        // support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
        for(int sy =std::max(0, iap[1]-1); sy <std::min(iap[1]+3, BS[1]); ++sy)
        for(int sx =std::max(0, iap[0]-1); sx <std::min(iap[0]+3, BS[0]); ++sx)
        {
          Real p[3]; //info.pos(p, sx-1, sy-1, 0-1);
          p[0] = ox + h * (sx - 1 + 0.5);
          p[1] = oy + h * (sy - 1 + 0.5);
          p[2] = oz + h * (0  - 1 + 0.5);

          const Real dist0 = eulerDistSq2D(p, myP);

          changeFromComputationalFrame(p);
          #ifndef NDEBUG // check that change of ref frame does not affect dist
            const Real p0[3] = {rX[ss] +width[ss]*signp*norX[ss],
                                rY[ss] +width[ss]*signp*norY[ss], 0
            };
            const Real distC = eulerDistSq2D(p, p0);
            assert(std::fabs(distC-dist0)<2.2e-16);
          #endif
          const Real distP = eulerDistSq2D(p,pP), distM = eulerDistSq2D(p,pM);

          int close_s = ss, secnd_s = ss + (distP<distM? 1 : -1);
          Real dist1 = dist0, dist2 = distP<distM? distP : distM;
          if(distP < dist0 || distM < dist0) { // switch nearest surf point
            dist1 = dist2; dist2 = dist0;
            close_s = secnd_s; secnd_s = ss;
          }

          const Real dSsq = std::pow(rX[close_s]-rX[secnd_s], 2)
                           +std::pow(rY[close_s]-rY[secnd_s], 2);
          assert(dSsq > 2.2e-16);
          const Real cnt2ML = std::pow( width[close_s],2);
          const Real nxt2ML = std::pow( width[secnd_s],2);

          Real sign2d = 0;
          if(dSsq>=std::fabs(cnt2ML-nxt2ML))
          { // if no abrupt changes in width we use nearest neighbour
            const Real xMidl[3] = {rX[close_s], rY[close_s], 0};
            const Real grd2ML = eulerDistSq2D(p, xMidl);
            sign2d = grd2ML > cnt2ML ? -1 : 1;
          } else {
            // else we model the span between ellipses as a spherical segment
            // http://mathworld.wolfram.com/SphericalSegment.html
            const Real corr = 2*std::sqrt(cnt2ML*nxt2ML);
            const Real Rsq = (cnt2ML +nxt2ML -corr +dSsq) //radius of the sphere
                            *(cnt2ML +nxt2ML +corr +dSsq)/4/dSsq;
            const Real maxAx = std::max(cnt2ML, nxt2ML);
            const int idAx1 = cnt2ML> nxt2ML? close_s : secnd_s;
            const int idAx2 = idAx1==close_s? secnd_s : close_s;
            // 'submerged' fraction of radius:
            const Real d = std::sqrt((Rsq - maxAx)/dSsq); // (divided by ds)
            // position of the centre of the sphere:
            const Real xMidl[3] = {rX[idAx1] +(rX[idAx1]-rX[idAx2])*d,
                                   rY[idAx1] +(rY[idAx1]-rY[idAx2])*d, 0};
            const Real grd2Core = eulerDistSq2D(p, xMidl);
            sign2d = grd2Core > Rsq ? -1 : 1; // as always, neg outside
          }

          //since naca extends over z axis, loop over all block
          for(int sz = 0; sz < BS[2]; ++sz) {
            const Real pZ = org[2] + h*sz;
            // positive inside negative outside ... as usual
            const Real distZ = height[ss] - std::fabs(position[2] - pZ);
            const Real signZ = (0 < distZ) - (distZ < 0);
            const Real dist3D = std::min(signZ*distZ*distZ, sign2d*dist1);

            if(std::fabs(SDFLAB[sz][sy][sx]) > dist3D) {
              SDFLAB[sz][sy][sx] = dist3D;
              const bool inRange =(sz - 1 >=0 && sz - 1 < ScalarBlock::sizeZ &&
                                   sy - 1 >=0 && sy - 1 < ScalarBlock::sizeY &&
                                   sx - 1 >=0 && sx - 1 < ScalarBlock::sizeX);
              if (inRange)
              {
                UDEF[sz-1][sy-1][sx-1][0] = udef[0];
                UDEF[sz-1][sy-1][sx-1][1] = udef[1];
                UDEF[sz-1][sy-1][sx-1][2] = udef[2];
                // not chi yet, just used for interpolating udef:
                CHI [sz-1][sy-1][sx-1] = 1;
              }
            }
          }
          // Not chi yet, I stored squared distance from analytical boundary
          // distSq is updated only if curr value is smaller than the old one
        }
      }
    }
  }
}

void PutNacaOnBlocks::constructInternl(const Real h, const Real ox, const Real oy, const Real oz, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB*>&vSegments) const
{
  Real org[3] = {ox-h,oy-h,oz-h};
  const Real invh = 1.0/h;
  const Real EPS = 1e-15;
  CHIMAT & __restrict__ CHI = defblock->chi;
  auto & __restrict__ SDFLAB = defblock->sdfLab;
  UDEFMAT & __restrict__ UDEF = defblock->udef;
  static constexpr int BS[3] = {ScalarBlock::sizeX+2, ScalarBlock::sizeY+2, ScalarBlock::sizeZ+2};

  // construct the deformation velocities (P2M with hat function as kernel)
  for(size_t i=0; i < vSegments.size(); ++i)
  {
  const int firstSegm = std::max(vSegments[i]->s_range.first,            1);
  const int lastSegm =  std::min(vSegments[i]->s_range.second, cfish->Nm-2);
  for(int ss=firstSegm; ss<=lastSegm; ++ss)
  {
    // P2M udef of a slice at this s
    const Real myWidth = cfish->width[ss], myHeight = cfish->height[ss];
    assert(myWidth > 0 && myHeight > 0);
    //here we process also all inner points. Nw to the left and right of midl
    // add xtension here to make sure we have it in each direction:
    const int Nw = std::floor(myWidth/h); //floor bcz we already did interior
    for(int iw = -Nw+1; iw < Nw; ++iw)
    {
      const Real offsetW = iw * h;
      Real xp[3] = { cfish->rX[ss] + offsetW*cfish->norX[ss],
                     cfish->rY[ss] + offsetW*cfish->norY[ss], 0
      };
      changeToComputationalFrame(xp);
      xp[0] = (xp[0]-org[0])*invh; // how many grid points from this block
      xp[1] = (xp[1]-org[1])*invh; // origin is this fishpoint located at?
      Real udef[3] = { cfish->vX[ss] + offsetW*cfish->vNorX[ss],
                       cfish->vY[ss] + offsetW*cfish->vNorY[ss], 0
      };
      changeVelocityToComputationalFrame(udef);
      const Real ap[2] = { std::floor(xp[0]), std::floor(xp[1]) };
      const int iap[2] = { (int)ap[0], (int)ap[1] };
      Real wghts[2][2]; // P2M weights
      for(int c=0; c<2; ++c) {
        const Real t[2] = { // we floored, hat between xp and grid point +-1
            std::fabs(xp[c] -ap[c]), std::fabs(xp[c] -(ap[c] +1))
        };
        wghts[c][0] = 1.0 - t[0];
        wghts[c][1] = 1.0 - t[1];
      }

      for(int idz=0; idz<BS[2]; ++idz)
      {
        const Real pZ = org[2] + h*idz;
        // positive inside negative outside ... as usual
        const Real distZ = myHeight - std::fabs(position[2] - pZ);
        static constexpr Real one = 1;
        const Real wz = .5 + std::min(one, std::max(distZ*invh, -one))/2;
        const Real signZ = (0 < distZ) - (distZ < 0);
        const Real distZsq = signZ*distZ*distZ;

        using std::max; using std::min;
        for(int sy=max(0,0-iap[1]); sy<min(2,BS[1]-iap[1]); ++sy)
        for(int sx=max(0,0-iap[0]); sx<min(2,BS[0]-iap[0]); ++sx) {
          const Real wxwywz = wz * wghts[1][sy] * wghts[0][sx];
          const int idx = iap[0]+sx, idy = iap[1]+sy;
          assert(idx>=0 && idx<BS[0]);
          assert(idy>=0 && idy<BS[1]);
          assert(wxwywz>=0 && wxwywz<=1);
          if (idz -1 >=0 && idz-1 < ScalarBlock::sizeZ &&
              idy -1 >=0 && idy-1 < ScalarBlock::sizeY &&
              idx -1 >=0 && idx-1 < ScalarBlock::sizeX)
          {
            UDEF[idz-1][idy-1][idx-1][0] += wxwywz*udef[0];
            UDEF[idz-1][idy-1][idx-1][1] += wxwywz*udef[1];
            UDEF[idz-1][idy-1][idx-1][2] += wxwywz*udef[2];
            CHI [idz-1][idy-1][idx-1] += wxwywz;
          }
          // set sign for all interior points:
          if(std::fabs(SDFLAB[idz][idy][idx]+1)<EPS) SDFLAB[idz][idy][idx] = distZsq;
        }
      }
    }
  }
  }
}

CubismUP_3D_NAMESPACE_END
