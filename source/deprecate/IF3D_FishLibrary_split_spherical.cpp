//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//


#include "IF3D_FishLibrary.h"

void FishMidlineData::writeMidline2File(const int step_id, std::string filename)
{
  char buf[500];
  sprintf(buf, "%s_midline_%07d.txt", filename, step_id);
  FILE * f = fopen(buf, "w");
  fprintf(f, "s x y vX vY\n");
  for (int i=0; i<Nm; i++) {
    //dummy.changeToComputationalFrame(temp);
    //dummy.changeVelocityToComputationalFrame(udef);
    fprintf(f, "%g %g %g %g %g\n", rS[i],rX[i],rY[i],vX[i],vY[i]);
  }
}

void FishMidlineData::_computeMidlineNormals()
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

Real FishMidlineData::integrateLinearMomentum(double CoM[2], double vCoM[2])
{   // already worked out the integrals for r, theta on paper
  // remaining integral done with composite trapezoidal rule
  // minimize rhs evaluations --> do first and last point separately
  double _vol=0, _cmx=0, _cmy=0, _lmx=0, _lmy=0;
  #pragma omp parallel for reduction(+:_vol,_cmx,_cmy,_lmx,_lmy)
  for(int i=0;i<Nm;++i) {
    const double ds = (i==0) ? rS[1]-rS[0] :
        ((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
    const double fac1 = _integrationFac1(i);
    const double fac2 = _integrationFac2(i);
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

  assert(vol> std::numeric_limits<double>::epsilon());
  const double ivol = 1.0/vol;

  CoM[0]*=ivol;
  CoM[1]*=ivol;
  vCoM[0]=linMom[0]*ivol;
  vCoM[1]=linMom[1]*ivol;
  //printf("%f %f %f %f %f\n",CoM[0],CoM[1],vCoM[0],vCoM[1], vol);
  return vol;
}

void FishMidlineData::integrateAngularMomentum(double& angVel)
{
  // assume we have already translated CoM and vCoM to nullify linear momentum
  // already worked out the integrals for r, theta on paper
  // remaining integral done with composite trapezoidal rule
  // minimize rhs evaluations --> do first and last point separately
  double _J = 0, _am = 0;
  #pragma omp parallel for reduction(+:_J,_am)
  for(int i=0;i<Nm;++i) {
    const double ds = (i==0) ? rS[1]-rS[0] :
        ((i==Nm-1) ? rS[Nm-1]-rS[Nm-2] :rS[i+1]-rS[i-1]);
    const double fac1 = _integrationFac1(i);
    const double fac2 = _integrationFac2(i);
    const double fac3 = _integrationFac3(i);
    const double tmp_M = (rX[i]*vY[i] - rY[i]*vX[i])*fac1
      + (rX[i]*vNorY[i] -rY[i]*vNorX[i] +vY[i]*norX[i] -vX[i]*norY[i])*fac2
      + (norX[i]*vNorY[i] - norY[i]*vNorX[i])*fac3;

    const double tmp_J = (rX[i]*rX[i] + rY[i]*rY[i])*fac1
      + 2*(rX[i]*norX[i] + rY[i]*norY[i])*fac2
      + fac3;

    _am += 0.5*tmp_M*ds;
    _J += 0.5*tmp_J*ds;
  }

  J=_J*M_PI;
  angMom=_am*M_PI;
  assert(J>std::numeric_limits<double>::epsilon());
  angVel = angMom/J;
}

void FishMidlineData::changeToCoMFrameLinear(const double CoM_internal[2], const double vCoM_internal[2])
{
  for(int i=0;i<Nm;++i) {
    rX[i]-=CoM_internal[0];
    rY[i]-=CoM_internal[1];
    vX[i]-=vCoM_internal[0];
    vY[i]-=vCoM_internal[1];
  }
}

void FishMidlineData::changeToCoMFrameAngular(const double theta_internal, const double angvel_internal)
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

void FishMidlineData::computeSurface()
{
  const int Nskin = lowerSkin->Npoints;
  // Compute surface points by adding width to the midline points
  #pragma omp parallel for
  for(int i=0; i<Nskin; ++i)
  {
    double norm[2] = {norX[iFishStart + i], norY[iFishStart + i]};
    double const norm_mod1 = std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] /= norm_mod1;
    norm[1] /= norm_mod1;
    assert(width[iFishStart + i] >= 0);
    lowerSkin->xSurf[i] = rX[iFishStart+i] - width[iFishStart+i]*norm[0];
    lowerSkin->ySurf[i] = rY[iFishStart+i] - width[iFishStart+i]*norm[1];
    upperSkin->xSurf[i] = rX[iFishStart+i] + width[iFishStart+i]*norm[0];
    upperSkin->ySurf[i] = rY[iFishStart+i] + width[iFishStart+i]*norm[1];
  }
}

void FishMidlineData::computeSkinNormals(const double theta_comp, const double CoM_comp[3])
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

    const double normL = std::sqrt( std::pow(lowerSkin->normXSurf[i],2) +
                                    std::pow(lowerSkin->normYSurf[i],2) );
    const double normU = std::sqrt( std::pow(upperSkin->normXSurf[i],2) +
                                    std::pow(upperSkin->normYSurf[i],2) );

    lowerSkin->normXSurf[i] /= normL;
    upperSkin->normXSurf[i] /= normU;
    lowerSkin->normYSurf[i] /= normL;
    upperSkin->normYSurf[i] /= normU;

    //if too close to the head or tail, consider a point further in, so that we are pointing out for sure
    const int ii = (i<8) ? 8 : ((i > Nskin-9) ? Nskin-9 : i);

    const Real dirL =
      lowerSkin->normXSurf[i] * (lowerSkin->midX[i]-rX[iFishStart+ii]) +
      lowerSkin->normYSurf[i] * (lowerSkin->midY[i]-rY[iFishStart+ii]);
    const Real dirU =
      upperSkin->normXSurf[i] * (upperSkin->midX[i]-rX[iFishStart+ii]) +
      upperSkin->normYSurf[i] * (upperSkin->midY[i]-rY[iFishStart+ii]);

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

void FishMidlineData::surfaceToCOMFrame(const double theta_internal, const double CoM_internal[2])
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
  }
}

void FishMidlineData::surfaceToComputationalFrame(const double theta_comp, const double CoM_interpolated[3])
{
  _prepareRotation2D(theta_comp);

  for(int i=0; i<upperSkin->Npoints; ++i)
  {
    _rotate2D(upperSkin->xSurf[i], upperSkin->ySurf[i]);
    upperSkin->xSurf[i] += CoM_interpolated[0];
    upperSkin->ySurf[i] += CoM_interpolated[1];
    _rotate2D(lowerSkin->xSurf[i], lowerSkin->ySurf[i]);
    lowerSkin->xSurf[i] += CoM_interpolated[0];
    lowerSkin->ySurf[i] += CoM_interpolated[1];
  }
}

void VolumeSegment_OBB::prepare(std::pair<int, int> _s_range, const Real bbox[3][2], const Real safe_dist)
{
  safe_distance = safe_dist;
  s_range.first = _s_range.first;
  s_range.second = _s_range.second;
  for(int i=0; i<3; ++i) {
    w[i] = 0.5*(bbox[i][1]-bbox[i][0]) + safe_distance;
    c[i] = 0.5*(bbox[i][1]+bbox[i][0]);
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

void VolumeSegment_OBB::changeToComputationalFrame(const double position[3], const double quaternion[4])
{
  // we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
  const Real a = quaternion[0];
  const Real x = quaternion[1];
  const Real y = quaternion[2];
  const Real z = quaternion[3];
  const Real Rmatrix[3][3] = {
      {1.-2*(y*y+z*z),    2*(x*y-z*a),    2*(x*z+y*a)},
      {   2*(x*y+z*a), 1.-2*(x*x+z*z),    2*(y*z-x*a)},
      {   2*(x*z-y*a),    2*(y*z+x*a), 1.-2*(x*x+y*y)}
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
      0.5*(end[0] - start[0]) + safe_distance,
      0.5*(end[1] - start[1]) + safe_distance,
      0.5*(end[2] - start[2]) + safe_distance
  };

  const Real AABB_c[3] = { //block center
    0.5*(end[0] + start[0]),
    0.5*(end[1] + start[1]),
    0.5*(end[2] + start[2])
  };

  const Real AABB_box[3][2] = {
    {AABB_c[0] - AABB_w[0],  AABB_c[0] + AABB_w[0]},
    {AABB_c[1] - AABB_w[1],  AABB_c[1] + AABB_w[1]},
    {AABB_c[2] - AABB_w[2],  AABB_c[2] + AABB_w[2]}
  };

  assert(AABB_w[0]>0 && AABB_w[1]>0 && AABB_w[2]>0);

  // Now Identify the ones that do not intersect
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

//inline void operator()(const BlockInfo& info, FluidBlock3D& b) const
void PutFishOnBlocks::operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const oblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  {
    const int N = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
    for(int i=0; i<N; ++i) *(&(oblock->chi[0][0][0])+i) = -1;
  }

  constructShape(info, b, oblock, vSegments);
  //constructDefVel(info, b, oblock, vSegments);
  //signedDistanceSqrt(info, b, oblock, vSegments);
}

void PutFishOnBlocks::constructShape(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint;
  const int iFishEnd = cfish->iFishEnd, iFishStart = cfish->iFishStart;
  const Real* const rX = cfish->rX;
  const Real* const rY = cfish->rY;
  const Real* const rS = cfish->rS;
  const Real* const norX = cfish->norX;
  const Real* const norY = cfish->norY;
  const Real* const width = cfish->width;
  const Real* const height = cfish->height;
  static constexpr int BS[3] = {FluidBlock::sizeX, FluidBlock::sizeY, FluidBlock::sizeZ};
  // construct the shape (P2M with min(distance) as kernel) onto defblocks
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
    //iterate over segments contained in the vSegm intersecting this block:
    const int firstSegm = max(vSegments[i].s_range.first,  iFishStart+1);
    const int lastSegm =  min(vSegments[i].s_range.second,   iFishEnd-1);
    for(int ss=firstSegm; ss<=lastSegm; ++ss) {
      assert(height[ss]>0 && width[ss]>0);
      // fill chi by crating an ellipse around ss and finding all near neighs
      // assume width is major axis, else correction:
      const Real offset = height[ss] > width[ss] ? M_PI/2 : 0;
      const Real ell_a = std::max(height[ss], width[ss]);
      // max distance between two points is ell_a * sin(dtheta): set it to dx/2
      const Real dtheta_tgt = std::fabs(std::asin(h/(ell_a+h)/2));
      const int Ntheta = std::ceil(2*M_PI/dtheta_tgt);
      const Real dtheta = 2*M_PI/((Real) Ntheta);

      for(int tt=0; tt<Ntheta; ++tt) {
        const Real theta = tt*dtheta + offset;
        const Real sinth = std::sin(theta), costh = std::cos(theta);
        // create a surface point
        Real myP[3] ={rX[ss+0] +width[ss+0]*costh*norX[ss+0],
                      rY[ss+0] +width[ss+0]*costh*norY[ss+0], height[ss+0]*sinth
        };
        changeToComputationalFrame(myP);
        const Real pP[3] = {rX[ss+1] +width[ss+1]*costh*norX[ss+1],
                      rY[ss+1] +width[ss+1]*costh*norY[ss+1], height[ss+1]*sinth
        };
        const Real pM[3] = {rX[ss-1] +width[ss-1]*costh*norX[ss-1],
                      rY[ss-1] +width[ss-1]*costh*norY[ss-1], height[ss-1]*sinth
        };
        // myP is now lab frame, find index of the fluid elem near it
        const int iap[3] = {
            std::floor((myP[0]-org[0])*invh),
            std::floor((myP[1]-org[1])*invh),
            std::floor((myP[2]-org[2])*invh)
        };

        // support is two points left, two points right --> Towers Chi
        // will be one point left, one point right, but needs SDF wider
        for(int sz =std::max(0, iap[2]-1); sz <std::min(iap[2]+3, BS[2]); ++sz)
        for(int sy =std::max(0, iap[1]-1); sy <std::min(iap[1]+3, BS[1]); ++sy)
        for(int sx =std::max(0, iap[0]-1); sx <std::min(iap[0]+3, BS[0]); ++sx)
        {
          Real p[3];
          info.pos(p, sx, sy, sz);
          const Real dist0 = eulerDistSq3D(p, myP);
          assert(dist0 <= 12*h*h);
          if(std::fabs(defblock->chi[sz][sy][sx]) < dist0) continue;

          changeFromComputationalFrame(p);
          #ifndef NDEBUG // check that change of ref frame does not affect dist
            Real p0[3] = {rX[ss] +width[ss]*costh*norX[ss],
                          rY[ss] +width[ss]*costh*norY[ss], height[ss]*sinth
            };
            const Real distC = eulerDistSq3D(p, p0);
            assert(std::fabs(distC-dist0)<2.2e-16);
          #endif
          const Real distP = eulerDistSq3D(p,pP), distM = eulerDistSq3D(p,pM);

          int close_s = ss;
          const int dir = distP<distM? 1 : -1;
          Real distSq = dist0;
          if(distP < distSq || distM < distSq) { // find nearest surf point
            distSq = distP<distM? distP : distM;
            close_s = ss + dir;
            /*
            while (0) {
              const int test_s = close_s + dir;
              if(test_s>iFishEnd || test_s<iFishStart) break;
              const Real s[3] = { rX[test_s]+width[test_s]*costh*norX[test_s],
               rY[test_s]+width[test_s]*costh*norY[test_s], height[test_s]*sinth
              };
              const Real distN = eulerDistSq3D(p, s);
              if(distN < distSq) {
                close_s = test_s;
                distSq = distN;
              } else break;
            }
            */
          }
          assert(distSq <= dist0);
          assert(close_s <= iFishEnd && close_s >= iFishStart);

          defblock->sectionMarker[sz][sy][sx] = cfish->rS[close_s];
          if (close_s < iFishEnd && close_s > iFishStart) // normal execution
          { // perform spherical integration between two nearest surf points
            // compute the two distances for the 2 near neighs
            const Real cnt2ML = std::pow( width[ss]*costh,2)
                               +std::pow(height[ss]*sinth,2);
            const Real nxt2ML = std::pow( width[ss+dir]*costh,2)
                               +std::pow(height[ss+dir]*sinth,2);
            const Real distN = distP<distM? distP:distM, W = 1/(distN+dist0);
            // spherical interp: surf to midline and closest midline point
            const Real surf2ML = (cnt2ML*distN +nxt2ML*dist0)*W;
            const Real xMidl[3] ={ (rX[ss]*distN + rX[ss+dir]*dist0)*W,
                                   (rY[ss]*distN + rY[ss+dir]*dist0)*W, 0};
            const Real grid2ML = eulerDistSq3D(p, xMidl);
            const Real sign = grid2ML > surf2ML ? -1 : 1;
            defblock->chi[sz][sy][sx] = sign*distSq;
            continue;
          }

          //we treat tips as spherical caps (wikipedia to find sphere centre)
          assert(close_s == iFishEnd || close_s == iFishStart);
          // index of nearest inner point:
          const int ints = std::max(iFishStart+1, std::min(iFishEnd-1,close_s));
          const Real dSsq = std::pow(rX[ints]-rX[close_s], 2)
                           +std::pow(rY[ints]-rY[close_s], 2);
          assert(dSsq > 2.2e-16);
          const Real dWsq = std::pow( width[ints]*costh,2)
                           +std::pow(height[ints]*sinth,2);
          // formula of radius given edge (width/height) and height (ds) of cap
          const Real srf2Core = (dWsq + dSsq)/dSsq/2; // here divided by dS
          const Real dist2int = srf2Core - 1; // 'submerged' fraction of radius
          if(dist2int > 0) {
            const Real xMidl[3] = {rX[ints] +(rX[ints]-rX[close_s])*dist2int,
                                rY[ints] +(rY[ints]-rY[close_s])*dist2int, 0};
            const Real grd2Core = eulerDistSq3D(p, xMidl);
            cout<<dist2int<<" "<<xMidl[0]<<" "<<xMidl[1]<<" "<<dSsq<<" "<<grd2Core<<endl; fflush(0);
            const Real sign = grd2Core > srf2Core*srf2Core*dSsq ? -1 : 1;
            defblock->chi[sz][sy][sx] = sign*distSq;
          }
          else // else the tip is sharp (tail) rather than blunt (head)
          {    // spherical approx is wrong and we just use closest ellipse
            const Real xMidl[3] = {rX[ints], rY[ints], 0};
            const Real srf2ML = std::pow( width[ints]*costh,2)
                               +std::pow(height[ints]*sinth,2);
            const Real grd2ML = eulerDistSq3D(p, xMidl);
            const Real sign = grd2ML > srf2ML ? -1 : 1;
            defblock->chi[sz][sy][sx] = sign*distSq;
          }
          // Not chi yet, I stored squared distance from analytical boundary
          // distSq is updated only if curr value is smaller than the old one
        }
      }
    }
  }
}

void PutFishOnBlocks::constructDefVel(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint;
  // construct the deformation velocities (P2M with hat function as kernel)
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
  for(int ss=vSegments[i].s_range.first; ss<=vSegments[i].s_range.second; ++ss)
  {
    assert(ss>=0 && ss<=cfish->Nm-1);
    // P2M udef of a slice at this s
    const Real myWidth =(ss<cfish->iFishStart? cfish->width[cfish->iFishStart]
                       :(ss>cfish->iFishEnd  ? cfish->width[cfish->iFishEnd]
                       : cfish->width[ss]));
    const Real myHeight=(ss<cfish->iFishStart?cfish->height[cfish->iFishStart]
                       :(ss>cfish->iFishEnd  ?cfish->height[cfish->iFishEnd]
                       : cfish->height[ss]));
    // towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice
    const Real extension = NPPEXT*h; // (and so we can floor!) G changed to 2

    const int Nh = std::floor( (myHeight+extension)/h );
    for(int ih=-Nh; ih<=Nh; ++ih)
    {
      const Real offsetH = ih * h;
      // add an extra extension when width == 0
      // (to deal with large curvatures near head and/or tail):
      const Real currentWidth = myWidth==0 ? extension
           : myWidth * std::sqrt(1 -std::pow(offsetH/(myHeight+extension), 2));
      const Real actualWidth =
          cfish->height[ss]<=0 or std::fabs(offsetH)>=cfish->height[ss] ? 0
        : cfish->width[ss]* std::sqrt(1-std::pow(offsetH/cfish->height[ss], 2));
      // add xtension here to make sure we have it in each direction:
      const int Nw = std::floor( (currentWidth+extension)/h );
      for(int iw=-Nw; iw<=Nw; ++iw)
      {
        const Real offsetW = iw * h;
        Real xp[3] = { cfish->rX[ss] + offsetW*cfish->norX[ss],
            cfish->rY[ss] + offsetW*cfish->norY[ss], offsetH
        };
        changeToComputationalFrame(xp);
        xp[0] = (xp[0]-org[0])*invh; // how many grid points
        xp[1] = (xp[1]-org[1])*invh; // from this block origin
        xp[2] = (xp[2]-org[2])*invh; // is this fishpoint located at?
        Real udef[3] = { cfish->vX[ss] + offsetW*cfish->vNorX[ss],
            cfish->vY[ss] + offsetW*cfish->vNorY[ss], 0
        };
        changeVelocityToComputationalFrame(udef);
        const Real ap[3] = {
            std::floor(xp[0]), std::floor(xp[1]), std::floor(xp[2])
        };
        const int iap[3] = { (int)ap[0], (int)ap[1], (int)ap[2] };
        Real wghts[3][2]; // P2M weights
        for(int c=0; c<3; ++c) {
          const Real t[2] = { // we floored, hat between xp and grid point +-1
              std::fabs(xp[c] -ap[c]), std::fabs(xp[c] -(ap[c] +1))
          };
          wghts[c][0] = 1.0 - t[0];
          wghts[c][1] = 1.0 - t[1];
        }
        const bool isInside =  (std::fabs(offsetW) < actualWidth)
                            && (std::fabs(offsetH) < cfish->height[ss]);
        for(int sz=max(0,0-iap[2]); sz<min(2,FluidBlock::sizeZ-iap[2]); ++sz) {
        for(int sy=max(0,0-iap[1]); sy<min(2,FluidBlock::sizeY-iap[1]); ++sy) {
        for(int sx=max(0,0-iap[0]); sx<min(2,FluidBlock::sizeX-iap[0]); ++sx) {
          const Real wxwywz = wghts[2][sz] * wghts[1][sy] * wghts[0][sx];
          assert(wxwywz>=0 && wxwywz<=1);
          const int idx[3] = { iap[0] + sx, iap[1] + sy, iap[2] + sz };
          assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
          assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
          assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);
          defblock->udef[idx[2]][idx[1]][idx[0]][0] += wxwywz*udef[0];
          defblock->udef[idx[2]][idx[1]][idx[0]][1] += wxwywz*udef[1];
          defblock->udef[idx[2]][idx[1]][idx[0]][2] += wxwywz*udef[2];
          b(idx[0],idx[1],idx[2]).tmpU += wxwywz;
          // set sign for all interior points
          if( (std::fabs(defblock->chi[idx[2]][idx[1]][idx[0]] +1) <
                5*std::numeric_limits<Real>::epsilon()) && isInside)
          defblock->chi[idx[2]][idx[1]][idx[0]] = 1.0;
      }
      }
      }
    }
    }
    }
  }
}

void PutFishOnBlocks::signedDistanceSqrt(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  // finalize signed distance function in tmpU
  const Real eps = std::numeric_limits<Real>::epsilon();
  for(int iz=0; iz<FluidBlock::sizeZ; iz++)
  for(int iy=0; iy<FluidBlock::sizeY; iy++)
  for(int ix=0; ix<FluidBlock::sizeX; ix++) {
    const Real normfac = b(ix,iy,iz).tmpU > eps ? b(ix,iy,iz).tmpU : 1.;
    defblock->udef[iz][iy][ix][0] /= normfac;
    defblock->udef[iz][iy][ix][1] /= normfac;
    defblock->udef[iz][iy][ix][2] /= normfac;
    // change from signed squared distance function to normal sdf
    b(ix,iy,iz).tmpU = defblock->chi[iz][iy][ix] > (Real)0 ?
      sqrt( defblock->chi[iz][iy][ix]) : -sqrt(-defblock->chi[iz][iy][ix]);
    //b(ix,iy,iz).tmpV = defblock->udef[iz][iy][ix][0]; //for debug
    //b(ix,iy,iz).tmpW = defblock->udef[iz][iy][ix][1]; //for debug

    // All points that are not chi=0 in the targeted section, are captured here. When we loop through SurfaceBlocks for computing torque, the extraneous points captured here will be left out, so hakunamatata.
    defblock->sectionMarker[iz][iy][ix] *= std::fabs(defblock->chi[iz][iy][ix]) > 0 ? 1 : 0;
  }
}

void PutNacaOnBlocks::signedDistanceSqrt(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  // finalize signed distance function in tmpU
  const Real eps = std::numeric_limits<Real>::epsilon();
  for(int iz=0; iz<FluidBlock::sizeZ; iz++)
  for(int iy=0; iy<FluidBlock::sizeY; iy++)
  for(int ix=0; ix<FluidBlock::sizeX; ix++) {
    // change from signed squared distance function to normal sdf
    b(ix,iy,iz).tmpU = defblock->chi[iz][iy][ix] > (Real)0 ?
      sqrt( defblock->chi[iz][iy][ix]) : -sqrt(-defblock->chi[iz][iy][ix]);
    //b(ix,iy,iz).tmpV = defblock->udef[iz][iy][ix][0]; //for debug
    //b(ix,iy,iz).tmpW = defblock->udef[iz][iy][ix][1]; //for debug
  }
}

void PutNacaOnBlocks::constructShape(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real invh = 1.0/info.h_gridpoint;
  // construct the shape (P2M with min(distance) as kernel) onto defblocks
  for(int i=0;i<(int)vSegments.size();++i) {
    //iterate over segments contained in the vSegm intersecting this block:
    const int firstSegm = max(vSegments[i].s_range.first,  cfish->iFishStart);
    const int lastSegm  = min(vSegments[i].s_range.second, cfish->iFishEnd);
    for(int ss=firstSegm; ss<=lastSegm; ++ss) {
      assert(ss>=cfish->iFishStart && ss<=cfish->iFishEnd);
      const Real myP_distSqPlane = std::pow(cfish->width[ss], 2);
      //for each segment, we have one point to left and right of midl
      for(int signp = -1; signp <= 1; signp+=2) {
        // create a surface point
        // special treatment of tail (width = 0 --> no ellipse, just line)
        Real myP[3] = { cfish->rX[ss] +cfish->width[ss]*signp*cfish->norX[ss],
          cfish->rY[ss] +cfish->width[ss]*signp*cfish->norY[ss], 0
        };
        changeToComputationalFrame(myP);
        const int iap[2] = {
            (int)std::floor((myP[0]-org[0])*invh),
            (int)std::floor((myP[1]-org[1])*invh)
        };
        // support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
        const int start[2] = { std::max(-1, 0-iap[0]), std::max(-1, 0-iap[1]) };
        const int end[3] = {
            std::min(+3, FluidBlock::sizeX - iap[0]),
            std::min(+3, FluidBlock::sizeY - iap[1])
        };

        for(int sy=start[1]; sy<end[1];++sy)
        for(int sx=start[0]; sx<end[0];++sx) {
          const int idx = iap[0] + sx, idy = iap[1] + sy;
          assert(idx>=0&&idx<FluidBlock::sizeX&&idy>=0&&idy<FluidBlock::sizeY);
          //since naca extends over z axis, loop over all block
          for(int idz = 0; idz < FluidBlock::sizeZ; ++idz)
          {
            Real p[3];
            info.pos(p, idx, idy, idz);
            //if(std::fabs(p[2]-position[2])>cfish->height[ss]) continue;

            const Real diff[2] = {p[0]-myP[0], p[1]-myP[1]};
            const Real distSq  = std::pow(diff[0],2)+std::pow(diff[1],2);
            int closest;
            const Real distToMidlineSq=getSmallerDistToMidLPlanar(ss,p,closest);

            //if(myP_distSqPlane < std::numeric_limits<Real>::epsilon()) {
            //  const Real distH = std::sqrt(std::pow(p[0]-cfish->rX[closest],2)
            //                              +std::pow(p[1]-cfish->rY[closest],2));
            //  const Real sign = distH>cfish->width[closest] ? -1 : 1;
              // Not chi yet, here stored squared distance from analytical boundary
              // Update the distSquared value only if the computed value is smaller than the old one
            //  defblock->chi[idz][idy][idx] =
            //      (std::abs(defblock->chi[idz][idy][idx]) > distSq) ?
            //      sign*distSq : defblock->chi[idz][idy][idx];
            //} else {
              const Real sign = distToMidlineSq >= myP_distSqPlane ? -1 : 1;
              defblock->chi[idz][idy][idx] =
                  (std::abs(defblock->chi[idz][idy][idx]) > distSq) ?
                  sign*distSq : defblock->chi[idz][idy][idx];
            //}
          }
        }
      }
    }
  }
}

void PutNacaOnBlocks::constructDefVel(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  {
    const size_t N = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX*3;
    memset(&(defblock->udef[0][0][0]), 0, N*sizeof(Real));
  }

  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real invh = 1.0/info.h_gridpoint;
  // construct the deformation velocities (P2M with hat function as kernel)
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
  for(int ss=vSegments[i].s_range.first; ss<=vSegments[i].s_range.second; ++ss)
  {
    assert(ss>=0 && ss<=cfish->Nm-1);
    // P2M udef of a slice at this s
    const Real myWidth =(ss<cfish->iFishStart? cfish->width[cfish->iFishStart]
                       :(ss>cfish->iFishEnd  ? cfish->width[cfish->iFishEnd]
                       : cfish->width[ss]));
    //const Real myHeight=(ss<cfish->iFishStart?cfish->height[cfish->iFishStart]
    //                   :(ss>cfish->iFishEnd  ?cfish->height[cfish->iFishEnd]
    //                   : cfish->height[ss]));
    // towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice (and so we can floor!)
    const Real extension = NPPEXT*info.h_gridpoint; //G tmp changed back to 2
    // add an extra extension when width == 0
    // (to deal with large curvatures near head and/or tail):
    const Real currentWidth = myWidth<=0 ? extension : myWidth;
    //here we process also all inner points. Nw to the left and right of midl
    // add xtension here to make sure we have it in each direction:
    const int Nw = std::floor( (currentWidth+extension)/info.h_gridpoint);
    for(int iw = -Nw; iw <= Nw; ++iw)
    {
      const Real offsetW = iw*info.h_gridpoint;
      Real xp[3] = { cfish->rX[ss] + offsetW*cfish->norX[ss],
          cfish->rY[ss] + offsetW*cfish->norY[ss], 0
      };
      changeToComputationalFrame(xp);
      xp[0] = (xp[0]-org[0])*invh;
      xp[1] = (xp[1]-org[1])*invh;
      const Real ap[2] = { std::floor((Real)xp[0]), std::floor((Real)xp[1]) };
      const int iap[2] = { (int)ap[0], (int)ap[1] };
      // now we P2M
      const int start[2] = { std::max(0,0-iap[0]), std::max(0,0-iap[1]) };
      const int end[2] = {
          std::min(+2, FluidBlock::sizeX - iap[0]),
          std::min(+2, FluidBlock::sizeY - iap[1])
      };

      for(int sy=start[1]; sy<end[1];++sy) {
      for(int sx=start[0]; sx<end[0];++sx) {
        const int idx[3] = {
            iap[0] + sx,
            iap[1] + sy
        };
        assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
        assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);

        for(int sz = 0; sz < FluidBlock::sizeZ; ++sz)
        {
          const bool isInside = std::fabs(offsetW)+info.h_gridpoint/2 < myWidth;
          // set sign for all interior points
          if( std::fabs(defblock->chi[sz][idx[1]][idx[0]] + 1) <
                std::numeric_limits<Real>::epsilon() && isInside )
          defblock->chi[sz][idx[1]][idx[0]] = 1.0;
        }
      }
      }
    }
  }
  }
}

Real PutNacaOnBlocks::getSmallerDistToMidLPlanar(const int start_s, const Real x[3], int & final_s) const
{
  Real relX[3] = {x[0],x[1],x[2]};
  changeFromComputationalFrame(relX);

  const Real curDistSq =  std::pow(relX[0]-cfish->rX[start_s],2)
  + std::pow(relX[1]-cfish->rY[start_s],2);

  Real distSq;

  distSq = curDistSq; // check right
  const int sRight = find_closest_dist_planar(start_s, +1, relX, distSq);

  distSq = curDistSq; // check left
  const int sLeft = find_closest_dist_planar(start_s, -1, relX, distSq);

  if(sRight==start_s and sLeft==start_s) {
    final_s = start_s;
    return distSq;
  }

  // With a hinged tail, the following overzealous assert will catch, outputted problematique and confirmed in Matlab
  //assert(sRight==start_s or sLeft==start_s);

  int curr_s = start_s;
  int new_s = sRight == start_s ? sLeft : sRight;
  const int dir = new_s-curr_s;
  while(curr_s not_eq new_s) {
    curr_s = new_s;
    new_s = find_closest_dist_planar(curr_s, dir, relX, distSq);
  }

  final_s = new_s;
  return distSq;
}

void MidlineShapes::integrateBSpline(const double*const xc,
  const double*const yc, const int n, const double length,
  Real*const rS, Real*const res, const int Nm)
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
  for(int i=0; i<Nm; ++i) {
    res[i] = 0;
    if (rS[i]>0 and rS[i]<length) {
      const double dtt = 0.1*(rS[i]-rS[i-1]);
      while (true) {
        double xi = 0;
        gsl_bspline_eval(ti, B, bw);
        for (int j=0; j<n; j++) xi += xc[j]*gsl_vector_get(B, j);
        if (xi >= rS[i]) break;
        ti += dtt;
      }

      for (int j=0; j<n; j++) res[i] += yc[j]*gsl_vector_get(B, j);
    }
  }
  gsl_bspline_free(bw);
  gsl_vector_free(B);
}

void MidlineShapes::naca_width(const double t_ratio, const double L,
  Real*const rS, Real*const res, const int Nm)
{
  const Real a =  0.2969;
  const Real b = -0.1260;
  const Real c = -0.3516;
  const Real d =  0.2843;
  const Real e = -0.1015;
  const Real t = t_ratio*L;

  for(int i=0; i<Nm; ++i)
  {
    if ( rS[i]<=0 or rS[i]>=L ) res[i] = 0;
    else {
      const Real p = rS[i]/L;
      res[i] = 5*t* (a*std::sqrt(p) +b*p +c*p*p +d*p*p*p + e*p*p*p*p);
      /*
      if(s>0.99*L){ // Go linear, otherwise trailing edge is not closed - NACA analytical's fault
        const Real temp = 0.99;
        const Real y1 = 5*t* (a*sqrt(temp) +b*temp +c*temp*temp +d*temp*temp*temp + e*temp*temp*temp*temp);
        const Real dydx = (0-y1)/(L-0.99*L);
        return y1 + dydx * (s - 0.99*L);
      }else{ // NACA analytical
        return 5*t* (a*sqrt(p) +b*p +c*p*p +d*p*p*p + e*p*p*p*p);
      }
      */
    }
  }
}

void MidlineShapes::computeWidthsHeights(const string heightName,
  const string widthName, const double L, Real* const rS,
  Real* const height, Real* const width, const int nM, const int mpirank)
{
  if ( heightName.compare("largefin") == 0 ) {
    if(!mpirank)
      cout<<"Building object's height according to 'largefin' profile."<<endl;
    double xh[8] = {0, 0, .2*L, .4*L, .6*L, .8*L, L, L};
    double yh[8] = {0, .055*L, .18*L, .2*L, .064*L, .002*L, .325*L, 0};
    // TODO read second to last number from factory
    integrateBSpline(xh, yh, 8, L, rS, height, nM);
  } else
  if ( heightName.compare("tunaclone") == 0 ) {
    if(!mpirank)
      cout<<"Building object's height according to 'tunaclone' profile."<<endl;
    double xh[9] = {0, 0, 0.2*L, .4*L, .6*L, .9*L, .96*L, L, L};
    double yh[9] = {0, .05*L, .14*L, .15*L, .11*L, 0, .2*L, .23*L, 0};
    // i also found 0.14 instead of 0.2 for yh[6]. Or ... :
    //const Real yh[9] = {0, 5e-2*length, 1.4e-1*length, 1.5e-1*length, //1.1e-1*length, .0*length, 0.1*length, 0.2*length, 0};
    integrateBSpline(xh, yh, 9, L, rS, height, nM);
  } else
  if ( heightName.compare(0, 4, "naca") == 0 ) {
    double t_naca = std::stoi( string(heightName, 5), nullptr, 10 ) * 0.01;
    if(!mpirank)
      cout<<"Building object's height according to naca profile with adim. thickness param set to "<<t_naca<<" ."<<endl;
    naca_width(t_naca, L, rS, height, nM);
  } else {
    if(!mpirank)
      cout<<"Building object's height according to baseline profile."<<endl;
    double xh[8] = {0, 0, .2*L, .4*L, .6*L, .8*L, L, L};
    double yh[8] = {0, .055*L, .068*L, .076*L, .064*L, .0072*L, .11*L, 0};
    integrateBSpline(xh, yh, 8, L, rS, height, nM);
  }

  if ( widthName.compare("fatter") == 0 ) {
    if(!mpirank)
      cout<<"Building object's width according to 'fatter' profile."<<endl;
    double xw[6] = {0, 0, L/3., 2*L/3., L, L};
    double yw[6] = {0, 8.9e-2*L, 7.0e-2*L, 3.0e-2*L, 2.0e-2*L, 0};
    integrateBSpline(xw, yw, 6, L, rS, width, nM);
  } else
  if ( widthName.compare(0, 4, "naca") == 0 ) {
    double t_naca = std::stoi( string(widthName, 5), nullptr, 10 ) * 0.01;
    if(!mpirank)
      cout<<"Building object's width according to naca profile with adim. thickness param set to "<<t_naca<<" ."<<endl;
    naca_width(t_naca, L, rS, width, nM);
  } else {
    if(!mpirank)
      cout<<"Building object's width according to baseline profile."<<endl;
    double xw[6] = {0, 0, L/3., 2*L/3., L, L};
    double yw[6] = {0, 8.9e-2*L, 1.7e-2*L, 1.6e-2*L, 1.3e-2*L, 0};
    integrateBSpline(xw, yw, 6, L, rS, width, nM);
  }

  if(!mpirank) {
    FILE * heightWidth;
    heightWidth = fopen("widthHeight.txt","w");
    for(int i=0; i<nM; ++i)
      fprintf(heightWidth,"%f \t %f \t %f \n", rS[i], width[i], height[i]);
    fclose(heightWidth);
  }
}

#if 0
inline Real _width(const Real s, const Real L)
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

inline Real _height(const Real s, const Real L)
{
  if(s<0 or s>L) return 0;
  const double a=0.51*L;
  const double b=0.08*L;
  return b*std::sqrt(1 - std::pow((s-a)/a,2));
}

class CarlingFishMidlineData : public FishMidlineData
{
 protected:
  //burst-coast:
  const Real tStart;
  Real t0, t1, t2, t3, lowestAmp=1e12;
  const Real fac, inv;
  //const Real sHinge, AhingeTheta, ThingeTheta, hingePhi;
  const Real sHinge, ThingeTheta;
  Real sLeft, sRight;
  Real AhingeTheta, hingePhi;
  const bool bBurst, bHinge;
  Real kSpring=0.0;
  const Real kMaxSpring=100.0; // Expect torque values on the order of 1e-5 at steady, and 1e-3 at startup
  Real thetaOld = 0.0, avgTorque = 0.0, runningTorque = 0.0, timeNminus = 0.0;
  int prevTransition = 0;

  Real aParabola, bParabola, cParabola;
  const bool quadraticAmplitude;
  const Real quadraticFactor = 0.1;

  inline Real rampFactorSine(const Real t, const Real T) const
  {
    return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
  }

  inline Real rampFactorVelSine(const Real t, const Real T) const
  {
    return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
  }

  inline Real getQuadAmp(const Real s, const Real L)
  {
    return s*s*quadraticFactor/L;
  }

  inline Real getArg(const Real s, const Real L, const Real t, const Real T, const Real phaseShift)
  {
    return 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
  }

  inline void computeParabolaParams(const Real yLeft, const Real yPrimeLeft, const Real yPrimeRight)
  {
    aParabola = (yPrimeLeft-yPrimeRight)/(2*(sLeft-sRight));
    bParabola = yPrimeRight - 2*aParabola*sRight;
    cParabola = yLeft - aParabola*sLeft*sLeft - bParabola*sLeft;
  }

  Real getJointParabola(const Real s, const Real L)
  {
    return aParabola*s*s + bParabola*s + cParabola;
  }

  Real midline(const Real s, const Real t, const Real L, const Real T,
    const Real phaseShift)// const
    //const Real phaseShift) const
  {
    //const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
    const Real arg = getArg(s, L, t, T, phaseShift);

    double yCurrent;
    if(quadraticAmplitude){
      //yCurrent = (s*s*quadraticFactor/L) *std::sin(arg);
      yCurrent = getQuadAmp(s,L) *std::sin(arg);
    } else {
      yCurrent = fac * (s + inv*L)*std::sin(arg);
    }

    // Just made the joint a whole lot more complicated. Now a smooth parabolic joint instead of a prick
    if(bHinge){
      if(not quadraticAmplitude) abort();
      if(s>=sLeft){
        const double yLeft = getQuadAmp(sLeft,L) * std::sin(getArg(sLeft,L,t,T,phaseShift));
        const double yPrimeLeft =
          (2*quadraticFactor*sLeft/L) * std::sin(getArg(sLeft,L,t,T,phaseShift))
         +getQuadAmp(sLeft,L)*std::cos(getArg(sLeft,L,t,T,phaseShift))*2*M_PI/(L*waveLength);

        const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
        const double yPrimeRight = std::sin(currentTheta);

        computeParabolaParams(yLeft,yPrimeLeft, yPrimeRight);

        yCurrent = getJointParabola(s,L);

        if(s>=sRight){
          const Real yRight = getJointParabola(sRight,L);
          yCurrent = yRight + yPrimeRight*(s-sRight);
        }
      }
      /*if(s>sHinge){
        double yNot;
        if(quadraticAmplitude){
          yNot =  (sHinge*sHinge*quadraticFactor/L)*std::sin(2.0*M_PI*(sHinge/(waveLength*L) - t/T + phaseShift));
        }else{
          yNot =  fac *  (sHinge + inv*L)*std::sin(2.0*M_PI*(sHinge/(waveLength*L) - t/T + phaseShift));
        }
        const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
        const double dydsNot = std::sin(currentTheta);
        yCurrent = yNot + dydsNot*(s-sHinge);
      }*/
    }
    return yCurrent;
  }

  inline Real midlineVel(const Real s, const Real t, const Real L, const Real T,
    const Real phaseShift)
    //const Real phaseShift) const
  {
    const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
    double velCurrent;
    if(quadraticAmplitude){
      //velCurrent = - (s*s*quadraticFactor/L)*(2.0*M_PI/T)*std::cos(arg);
      velCurrent = (-2.0*M_PI/T)*getQuadAmp(s,L)*std::cos(getArg(s, L, t, T, phaseShift));
    }else{
      velCurrent = - fac*(s + inv*L)*(2.0*M_PI/T)*std::cos(arg);
    }

    if(bHinge){
      if(not quadraticAmplitude) abort();
      if(s>=sLeft){
        //const double yLeft = getQuadAmp(sLeft,L) * std::sin(getArg(sLeft,L,t,T,phaseShift));
        const double yLeftDot  = (-2.0*M_PI/T) * getQuadAmp(sLeft,L)*std::cos(getArg(sLeft,L,t,T,phaseShift));
        //const double yPrimeLeft = (2*quadraticFactor*sLeft/L) * sin(getArg(sLeft,L,t,T,phaseShift))
        //  + getQuadAmp(sLeft,L)*cos(getArg(sLeft,L,t,T,phaseShift))*2.0*M_PI/(L*waveLength);
        const double yPrimeLeftDot = (2*quadraticFactor*sLeft/L) * (-2*M_PI/T) * std::cos(getArg(sLeft,L,t,T,phaseShift))
          + getQuadAmp(sLeft,L)*(2*M_PI/T)*sin(getArg(sLeft,L,t,T,phaseShift))*2.0*M_PI/(L*waveLength);

        const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
        const double currentThetaDot = AhingeTheta * (2*M_PI/ThingeTheta)*std::cos(2.0*M_PI*(t/ThingeTheta + hingePhi));
        //const double yPrimeRight = std::sin(currentTheta);
              const double yPrimeRightDot = std::cos(currentTheta)*currentThetaDot;

        const double aDot = (yPrimeLeftDot - yPrimeRightDot)/(2*(sLeft-sRight));
        const double bDot = yPrimeRightDot - 2*sRight*aDot;
        const double cDot = yLeftDot - sLeft*sLeft*aDot - sLeft*bDot;
        velCurrent = aDot*s*s + bDot*s + cDot;

        if(s>=sRight){
          //const Real yRight = getJointParabola(sRight,L);
          //yCurrent = yRight + yPrimeRight*(s-sRight);
          const Real yRightDot = aDot*sRight*sRight + bDot*sRight + cDot;
          velCurrent = yRightDot + yPrimeRightDot*(s-sRight);
        }
      }
      /*
      if(s>sHinge) {
        //const double yNot =  4./33 *  (sHinge + 0.03125*L)*std::sin(2.0*M_PI*(sHinge/L - t/T + phaseShift));
        double velNot;
        double velNot;
        if(quadraticAmplitude){
          velNot =  -2.0*M_PI/T * (sHinge*sHinge*quadraticFactor/L)
            *std::cos(2.0*M_PI*(sHinge/(L*waveLength) - t/T + phaseShift));
        }else{
          velNot =  -2.0*M_PI/T * fac *  (sHinge + inv*L)*
          std::cos(2.0*M_PI*(sHinge/(L*waveLength) - t/T + phaseShift));
        }
        const double currentTheta = AhingeTheta * std::sin(2.0*M_PI*(t/ThingeTheta + hingePhi));
        const double currentThetaDot = AhingeTheta * 2.0*M_PI/ThingeTheta * std::cos(2.0*M_PI*(t/ThingeTheta + hingePhi));
        const double dydsNotDT = std::cos(currentTheta)*currentThetaDot;
        velCurrent = velNot + dydsNotDT*(s-sHinge);
      }
      */
    }
    return velCurrent;
  }

  inline Real midlineBeC(const Real s, const Real t, const Real L, const Real T,
    const Real phaseShift, const Real f) const
  {
    const Real arg = 2.0*M_PI*(s/(waveLength*L) - t/T + phaseShift);
    return f * fac * (s + inv*L)*std::sin(arg);
  }

  inline Real midlineVelBeC(const Real s, const Real t, const Real L,
    const Real T, const Real phaseShift, const Real f, const Real df) const
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
    //           : (shift+1)*0.5 + phaseShift;
    const Real phase = 0.0;

    if (time<tcoast) {
      printf("NCCUCDC.\n");
      abort();
    } else if (time<tfreeze) {
      const Real d = (time-tcoast)/(tfreeze-tcoast);
      const std::pair<double, double> retVal = cubicHermite(1.0, lowestAmp, d);
      f = retVal.first;
      df = retVal.second/(tfreeze-tcoast);
      //f = 1 - 3*d*d + 2*d*d*d;
    } else if (time<tburst) {
      //f = 0.0;
      f = lowestAmp;
      df = 0.0;
    } else if (time<tswim) {
      const Real d = (time-tburst)/(tswim-tburst);
      const std::pair<double, double> retVal = cubicHermite(lowestAmp, 1.0, d);
      f = retVal.first;
      df = retVal.second/(tswim-tburst);
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
        if(not(dx>0))  vX[i] = 0.0;
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
  : FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(false), tStart(1.0e12), bHinge(false), sHinge(0.0), AhingeTheta(0.0), ThingeTheta(0.0), hingePhi(0.0), quadraticAmplitude(true)
  /*: FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), fac(_fac), inv(0.03125), bBurst(false), tStart(1.0e12), bHinge(false), sHinge(0.0), AhingeTheta(0.0), ThingeTheta(0.0), hingePhi(0.0), quadraticAmplitude(false)
  {
  }*/
  {
    // Now, read the optimization params (Ahinge, _phiHinge, tail size) from the params file
    {
      Real _Ahinge, _phiHinge;
      ifstream reader("hingedParams.txt");
      if (reader.is_open()) {
        reader >> _Ahinge;
        reader >> _phiHinge;
        reader >> waveLength;
        //reader >> finSize;
        //printf("Read numbers = %f, %f, %f, %f\n", _Ahinge, _phiHinge, waveLength, finSize);
        printf("Read numbers = %f, %f, %f\n", _Ahinge, _phiHinge, waveLength);
        if(reader.eof()){
          cout << "Insufficient number of parameters provided for hingedFin" << endl; fflush(NULL); abort();
        }
        reader.close();
      } else {
        cout << "Could not open the correct 'params'.txt file" << endl; fflush(NULL);
        abort();
      }
    }

    // FinSize has now been updated with value read from text file. Recompute heights to over-write with updated values
    printf("Overwriting default tail-fin size for Plain Carling:\n");
    _computeWidthsHeights();

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
    FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext), tStart(1e12), fac(_fac), inv(0.03125), sHinge(_sHinge), ThingeTheta(Tperiod), AhingeTheta(M_PI*_Ahinge/180.0), hingePhi(_phiHinge/360.0), bBurst(false), bHinge(true), quadraticAmplitude(true)
    {
    }

  // Had to get rid of default value of _fac=0.1212121212121212. Was creating confusion for compiler with plain fish constructor
  CarlingFishMidlineData(const int Nm, const Real length, const Real Tperiod,
      const Real phaseShift, const Real dx_ext, const Real _sHinge,
      const Real _fac, const bool _equalHeight):
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
        reader >> waveLength;
        /*
        if(_equalHeight){
          finSize = 0.3;
        }else{
          reader >> finSize;
        }
        */
        //printf("Read numbers = %f, %f, %f, %f\n", _Ahinge, _phiHinge, waveLength, finSize);
        printf("Read numbers = %f, %f, %f\n", _Ahinge, _phiHinge, waveLength);
        /*reader >> sHinge2;
        reader >> kSpring;*/
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
      /*
      sHinge2 *= length;
      // UnRescaling: to avoid CMA trouble
      kSpring *= 1.0e-4;
      */
      sLeft  = sHinge - 0.02*length;
      sRight = sHinge + 0.02*length;
    }

    // FinSize has now been updated with value read from text file. Recompute heights to over-write with updated values
    printf("Overwriting default tail-fin size for Hinged Carling:\n");
    _computeWidthsHeights();

  }

  void computeMidline(const Real time, const Real dt) override
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

#endif
