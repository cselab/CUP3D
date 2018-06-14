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
  sprintf(buf, "%s_midline_%07d.txt", filename.c_str(), step_id);
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
    double norm[2] = {norX[i], norY[i]};
    double const norm_mod1 = std::sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0] /= norm_mod1;
    norm[1] /= norm_mod1;
    assert(width[i] >= 0);
    lowerSkin->xSurf[i] = rX[i] - width[i]*norm[0];
    lowerSkin->ySurf[i] = rY[i] - width[i]*norm[1];
    upperSkin->xSurf[i] = rX[i] + width[i]*norm[0];
    upperSkin->ySurf[i] = rY[i] + width[i]*norm[1];
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
      lowerSkin->normXSurf[i] * (lowerSkin->midX[i]-rX[ii]) +
      lowerSkin->normYSurf[i] * (lowerSkin->midY[i]-rY[ii]);
    const Real dirU =
      upperSkin->normXSurf[i] * (upperSkin->midX[i]-rX[ii]) +
      upperSkin->normYSurf[i] * (upperSkin->midY[i]-rY[ii]);

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

  #pragma omp parallel for
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

  #pragma omp parallel for
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

void VolumeSegment_OBB::changeToComputationalFrame(const double position[3], const double quaternion[4])
{
  // we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
  const Real a = quaternion[0];
  const Real x = quaternion[1];
  const Real y = quaternion[2];
  const Real z = quaternion[3];
  const double Rmatrix[3][3] = {
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
  //std::chrono::time_point<std::chrono::high_resolution_clock> t0, t1, t2, t3;
  //t0 = std::chrono::high_resolution_clock::now();
  constructSurface(info, b, oblock, vSegments);
  //t1 = std::chrono::high_resolution_clock::now();
  constructInternl(info, b, oblock, vSegments);
  //t2 = std::chrono::high_resolution_clock::now();
  signedDistanceSqrt(info, b, oblock, vSegments);
  //t3 = std::chrono::high_resolution_clock::now();
  //printf("%g %g %g\n",std::chrono::duration<double>(t1-t0).count(),
  //                    std::chrono::duration<double>(t2-t1).count(),
  //                    std::chrono::duration<double>(t3-t2).count());
}

void PutFishOnBlocks::constructSurface(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint;
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
  static constexpr int BS[3] = {FluidBlock::sizeX, FluidBlock::sizeY, FluidBlock::sizeZ};
  // construct the shape (P2M with min(distance) as kernel) onto defblocks
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
    //iterate over segments contained in the vSegm intersecting this block:
    const int firstSegm = max(vSegments[i].s_range.first,            1);
    const int lastSegm =  min(vSegments[i].s_range.second, cfish->Nm-2);
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
            (int)std::floor((myP[0]-org[0])*invh),
            (int)std::floor((myP[1]-org[1])*invh),
            (int)std::floor((myP[2]-org[2])*invh)
        };
        Real udef[3] = { vX[ss+0] +width[ss+0]*costh*vNorX[ss+0],
                         vY[ss+0] +width[ss+0]*costh*vNorY[ss+0], 0
        };
        changeVelocityToComputationalFrame(udef);

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

          int close_s = ss, secnd_s = ss + (distP<distM? 1 : -1);
          Real dist1 = dist0, dist2 = distP<distM? distP : distM;
          if(distP < dist0 || distM < dist0) { // switch nearest surf point
            dist1 = dist2; dist2 = dist0;
            close_s = secnd_s; secnd_s = ss;
          }

          const Real dSsq = std::pow(rX[close_s]-rX[secnd_s], 2)
                           +std::pow(rY[close_s]-rY[secnd_s], 2);
          assert(dSsq > 2.2e-16);
          const Real cnt2ML = std::pow( width[close_s]*costh,2)
                             +std::pow(height[close_s]*sinth,2);
          const Real nxt2ML = std::pow( width[secnd_s]*costh,2)
                             +std::pow(height[secnd_s]*sinth,2);

          const Real corr = 2*std::sqrt(cnt2ML*nxt2ML);
          defblock->sectionMarker[sz][sy][sx] = close_s;
          defblock->udef[sz][sy][sx][0] = udef[0];
          defblock->udef[sz][sy][sx][1] = udef[1];
          defblock->udef[sz][sy][sx][2] = udef[2];
          b(sx,sy,sz).tmpU = 1;
          if(dSsq >= cnt2ML+nxt2ML -corr) // if ds > delta radius
          { // if no abrupt changes in width we use nearest neighbour
            const Real xMidl[3] = {rX[close_s], rY[close_s], 0};
            const Real grd2ML = eulerDistSq3D(p, xMidl);
            const Real sign = grd2ML > cnt2ML ? -1 : 1;
            defblock->chi[sz][sy][sx] = sign*dist1;
          }
          else
          {
            // else we model the span between ellipses as a spherical segment
            // http://mathworld.wolfram.com/SphericalSegment.html
            const Real Rsq = (cnt2ML +nxt2ML -corr +dSsq) //radius of the spere
                            *(cnt2ML +nxt2ML +corr +dSsq)/4/dSsq;
            const Real maxAx = std::max(cnt2ML, nxt2ML);
            const int idAx1 = cnt2ML> nxt2ML? close_s : secnd_s;
            const int idAx2 = idAx1==close_s? secnd_s : close_s;
            // 'submerged' fraction of radius:
            const Real d = std::sqrt((Rsq - maxAx)/dSsq); // (divided by ds)
            // position of the centre of the sphere:
            const Real xMidl[3] = {rX[idAx1] +(rX[idAx1]-rX[idAx2])*d,
                                   rY[idAx1] +(rY[idAx1]-rY[idAx2])*d, 0};
            const Real grd2Core = eulerDistSq3D(p, xMidl);
            const Real sign = grd2Core > Rsq ? -1 : 1;
            defblock->chi[sz][sy][sx] = sign*dist1;
          }
          // Not chi yet, I stored squared distance from analytical boundary
          // distSq is updated only if curr value is smaller than the old one
        }
      }
    }
  }
}

void PutFishOnBlocks::constructInternl(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint, EPS = 1e-15;
  // construct the deformation velocities (P2M with hat function as kernel)
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
  const int firstSegm = max(vSegments[i].s_range.first,            1);
  const int lastSegm =  min(vSegments[i].s_range.second, cfish->Nm-2);
  for(int ss=firstSegm; ss<=lastSegm; ++ss)
  {
    // P2M udef of a slice at this s
    const Real myWidth = cfish->width[ss], myHeight = cfish->height[ss];
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
          assert(wghts[c][0]>=0 && wghts[c][1]>=0);
        }

        for(int sz=max(0,0-iap[2]); sz<min(2,FluidBlock::sizeZ-iap[2]); ++sz)
        for(int sy=max(0,0-iap[1]); sy<min(2,FluidBlock::sizeY-iap[1]); ++sy)
        for(int sx=max(0,0-iap[0]); sx<min(2,FluidBlock::sizeX-iap[0]); ++sx) {
          const Real wxwywz = wghts[2][sz] * wghts[1][sy] * wghts[0][sx];
          const int idx = iap[0]+sx, idy = iap[1]+sy, idz = iap[2]+sz;
          assert(idx>=0 && idx<FluidBlock::sizeX);
          assert(idy>=0 && idy<FluidBlock::sizeY);
          assert(idz>=0 && idz<FluidBlock::sizeZ);
          assert(wxwywz>=0 && wxwywz<=1);

          defblock->udef[idz][idy][idx][0] += wxwywz*udef[0];
          defblock->udef[idz][idy][idx][1] += wxwywz*udef[1];
          defblock->udef[idz][idy][idx][2] += wxwywz*udef[2];
          b(idx,idy,idz).tmpU += wxwywz;

          // set sign for all interior points
          if( std::fabs(defblock->chi[idz][idy][idx]+1)<EPS )
            defblock->chi[idz][idy][idx] = 1.0;
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
  }
}

void PutNacaOnBlocks::constructSurface(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint;
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
  static constexpr int BS[3] = {FluidBlock::sizeX, FluidBlock::sizeY, FluidBlock::sizeZ};

  // construct the shape (P2M with min(distance) as kernel) onto defblocks
  for(int i=0; i<(int)vSegments.size(); ++i) {
    //iterate over segments contained in the vSegm intersecting this block:
    const int firstSegm = max(vSegments[i].s_range.first,            1);
    const int lastSegm =  min(vSegments[i].s_range.second, cfish->Nm-2);
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
          Real p[3];
          info.pos(p, sx, sy, 0);
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
          for(int sz = 0; sz < FluidBlock::sizeZ; ++sz) {
            const Real pZ = org[2] + h*sz;
            // positive inside negative outside ... as usual
            const Real distZ = height[ss] - std::fabs(position[2] - pZ);
            const Real signZ = (0 < distZ) - (distZ < 0);
            const Real dist3D = std::min(signZ*distZ*distZ, sign2d*dist1);

            if(std::fabs(defblock->chi[sz][sy][sx]) > dist3D) {
              defblock->sectionMarker[sz][sy][sx] = close_s;
              defblock->udef[sz][sy][sx][0] = udef[0];
              defblock->udef[sz][sy][sx][1] = udef[1];
              defblock->udef[sz][sy][sx][2] = udef[2];
              defblock->chi [sz][sy][sx] = dist3D;
              b(sx,sy,sz).tmpU = 1;
            }
          }
          // Not chi yet, I stored squared distance from analytical boundary
          // distSq is updated only if curr value is smaller than the old one
        }
      }
    }
  }
}

void PutNacaOnBlocks::constructInternl(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
{
  Real org[3];
  info.pos(org, 0, 0, 0);
  const Real h = info.h_gridpoint, invh = 1.0/info.h_gridpoint, EPS = 1e-15;
  // construct the deformation velocities (P2M with hat function as kernel)
  for(int i=0; i<(int)vSegments.size(); ++i)
  {
  const int firstSegm = max(vSegments[i].s_range.first,            1);
  const int lastSegm =  min(vSegments[i].s_range.second, cfish->Nm-2);
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

      for(int idz=0; idz<FluidBlock::sizeZ; ++idz)
      {
        const Real pZ = org[2] + h*idz;
        // positive inside negative outside ... as usual
        const Real distZ = myHeight - std::fabs(position[2] - pZ);
        static constexpr Real one = 1;
        const Real wz = .5 + std::min(one, std::max(distZ*invh, -one))/2;
        const Real signZ = (0 < distZ) - (distZ < 0);
        const Real distZsq = signZ*distZ*distZ;

        for(int sy=max(0,0-iap[1]); sy<min(2,FluidBlock::sizeY-iap[1]); ++sy)
        for(int sx=max(0,0-iap[0]); sx<min(2,FluidBlock::sizeX-iap[0]); ++sx) {
          const Real wxwywz = wz * wghts[1][sy] * wghts[0][sx];
          const int idx = iap[0]+sx, idy = iap[1]+sy;
          assert(idx>=0 && idx<FluidBlock::sizeX);
          assert(idy>=0 && idy<FluidBlock::sizeY);
          assert(wxwywz>=0 && wxwywz<=1);

          defblock->udef[idz][idy][idx][0] += wxwywz*udef[0];
          defblock->udef[idz][idy][idx][1] += wxwywz*udef[1];
          defblock->udef[idz][idy][idx][2] += wxwywz*udef[2];
          b(idx,idy,idz).tmpU += wxwywz;

          // set sign for all interior points
          if( std::fabs(defblock->chi[idz][idy][idx]+1)<EPS )
            defblock->chi[idz][idy][idx] = distZsq;
        }
      }
    }
  }
  }
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
      const double dtt = (rS[i]-rS[i-1])/1e3;
      while (true) {
        double xi = 0;
        gsl_bspline_eval(ti, B, bw);
        for (int j=0; j<n; j++) xi += xc[j]*gsl_vector_get(B, j);
        if (xi >= rS[i]) break;
        if(ti + dtt > len)  break;
        else ti += dtt;
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

void MidlineShapes::stefan_width(const double L, Real*const rS, Real*const res, const int Nm)
{
  const double sb = .04*L;
  const double st = .95*L;
  const double wt = .01*L;
  const double wh = .04*L;

  for(int i=0; i<Nm; ++i)
  {
    if(rS[i]<=0 or rS[i]>=L) res[i] = 0;
    else {
      const Real s = rS[i];
      res[i] = (s<sb ? std::sqrt(2.0*wh*s-s*s) :
               (s<st ? wh-(wh-wt)*std::pow((s-sb)/(st-sb),2) :
               (wt * (L-s)/(L-st))));
    }
  }
}

void MidlineShapes::stefan_height(const double L, Real*const rS, Real*const res, const int Nm)
{
  const double a=0.51*L;
  const double b=0.08*L;

  for(int i=0; i<Nm; ++i)
  {
    if(rS[i]<=0 or rS[i]>=L) res[i] = 0;
    else {
      const Real s = rS[i];
      res[i] = b*std::sqrt(1 - std::pow((s-a)/a,2));
    }
  }
}

void MidlineShapes::danio_width(const double L, Real*const rS, Real*const res, const int Nm)
{
	const int nBreaksW = 11;
	const double breaksW[nBreaksW] = {0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0};
	const double coeffsW[nBreaksW-1][4] = {
		{0.0015713,       2.6439,            0,       -15410},
		{ 0.012865,       1.4882,      -231.15,        15598},
		{0.016476,      0.34647,       2.8156,      -39.328},
		{0.032323,      0.38294,      -1.9038,       0.7411},
		{0.046803,      0.19812,      -1.7926,       5.4876},
		{0.054176,    0.0042136,     -0.14638,     0.077447},
		{0.049783,    -0.045043,    -0.099907,     -0.12599},
		{ 0.03577,     -0.10012,      -0.1755,      0.62019},
		{0.013687,      -0.0959,      0.19662,      0.82341},
		{0.0065049,     0.018665,      0.56715,       -3.781}
	};

	for(int i=0; i<Nm; ++i)
	{
		if ( rS[i]<=0 or rS[i]>=L ) res[i] = 0;
		else {

			const double sNormalized = rS[i]/L;

			// Find current segment
			int currentSegW = 1;
			while(sNormalized>=breaksW[currentSegW]) currentSegW++;
			currentSegW--; // Shift back to the correct segment
			//if(rS[i]==L) currentSegW = nBreaksW-2; Not necessary - excluded by the if conditional


			// Declare pointer for easy access
			const double *paramsW = coeffsW[currentSegW];
			// Reconstruct cubic
			const double xxW = sNormalized - breaksW[currentSegW];
			res[i] = L*(paramsW[0] + paramsW[1]*xxW + paramsW[2]*pow(xxW,2) + paramsW[3]*pow(xxW,3));
		}
	}
}

void MidlineShapes::danio_height(const double L, Real*const rS, Real*const res, const int Nm)
{
	const int nBreaksH = 15;
	const double breaksH[nBreaksH] = {0, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.8, 0.85, 0.87, 0.9, 0.993, 0.996, 0.998, 1};
	const double coeffsH[nBreaksH-1][4] = {
		{0.0011746,        1.345,   2.2204e-14,      -578.62},
		{0.014046,       1.1715,      -17.359,        128.6},
		{0.041361,     0.40004,      -1.9268 ,      9.7029},
		{0.057759,      0.28013,     -0.47141,     -0.08102},
		{0.094281,     0.081843,     -0.52002,     -0.76511},
		{0.083728,     -0.21798,     -0.97909,       3.9699},
		{0.032727,     -0.13323,       1.4028,       2.5693},
		{0.036002,      0.22441,       2.1736,      -13.194},
		{0.051007,      0.34282,      0.19446,       16.642},
		{0.058075,      0.37057,        1.193,      -17.944},
		{0.069781,       0.3937,     -0.42196,      -29.388},
		{0.079107,     -0.44731,      -8.6211,  -1.8283e+05},
		{0.072751,      -5.4355,      -1654.1,  -2.9121e+05},
		{0.052934,      -15.546,      -3401.4,   5.6689e+05}
	};

	for(int i=0; i<Nm; ++i)
	{
		if ( rS[i]<=0 or rS[i]>=L ) res[i] = 0;
		else {

			const double sNormalized = rS[i]/L;

			// Find current segment
			int currentSegH = 1;
			while(sNormalized>=breaksH[currentSegH]) currentSegH++;
			currentSegH--; // Shift back to the correct segment
			//if(rS[i]==L) currentSegH = nBreaksH-2; Not necessary - excluded by the if conditional

			// Declare pointer for easy access
			const double *paramsH = coeffsH[currentSegH];
			// Reconstruct cubic
			const double xxH = sNormalized - breaksH[currentSegH];
			res[i] = L*(paramsH[0] + paramsH[1]*xxH + paramsH[2]*pow(xxH,2) + paramsH[3]*pow(xxH,3));
		}
	}
}

void MidlineShapes::computeWidthsHeights(const string heightName,
  const string widthName, const double L, Real* const rS,
  Real* const height, Real* const width, const int nM, const int mpirank)
{
  printf("height = %s, width=%s\n", heightName.c_str(), widthName.c_str()); fflush(NULL);
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
  } else
  if ( heightName.compare("danio") == 0 ) {
    if(!mpirank)
      cout<<"Building object's height according to Danio (zebrafish) profile from Maertens2017 (JFM)"<<endl;
    danio_height(L, rS, height, nM);
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
  } else
  if ( widthName.compare("danio") == 0 ) {
    if(!mpirank)
      cout<<"Building object's width according to Danio (zebrafish) profile from Maertens2017 (JFM)"<<endl;
    danio_width(L, rS, width, nM);
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
