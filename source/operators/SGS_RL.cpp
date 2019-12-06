//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Operator.h"
#include "SGS_RL.h"
#include "smarties.h"
#include "../spectralOperators/SpectralManip.h"
#include "../spectralOperators/HITtargetData.h"

#include <functional>
CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

struct ActionInterpolator
{
  const int NX, NY, NZ;
  static constexpr int NB = CUP_BLOCK_SIZE;

  std::vector<std::vector<std::vector<double>>> actions =
    std::vector<std::vector<std::vector<double>>> (
      NZ, std::vector<std::vector<double>>(NY, std::vector<double>(NX, 0) ) );

  ActionInterpolator(int _NX, int _NY, int _NZ) : NX(_NX), NY(_NY), NZ(_NZ)
  {
  }

  double operator()(const int bix, const int biy, const int biz,
    const int ix, const int iy, const int iz) const
  {
    // linear interpolation betwen element's block (bix, biy, biz) and second
    // nearest. figure out which from element's index (ix, iy, iz) in block:
    const int nbix = ix < NB/2 ? bix - 1 : bix + 1;
    const int nbiy = iy < NB/2 ? biy - 1 : biy + 1;
    const int nbiz = iz < NB/2 ? biz - 1 : biz + 1;
    // distance from second nearest block along its direction:
    const Real dist_nbix = ix < NB/2 ? NB/2 + ix + 0.5 : 3*NB/2 - ix - 0.5;
    const Real dist_nbiy = iy < NB/2 ? NB/2 + iy + 0.5 : 3*NB/2 - iy - 0.5;
    const Real dist_nbiz = iz < NB/2 ? NB/2 + iz + 0.5 : 3*NB/2 - iz - 0.5;
    // distance from block's center:
    const Real dist_bix = std::fabs(ix + 0.5 - NB/2);
    const Real dist_biy = std::fabs(iy + 0.5 - NB/2);
    const Real dist_biz = std::fabs(iz + 0.5 - NB/2);

    double weighted_sum_act = 0;
    double sum_acts_weights = 0;
    for(int z = 0; z < 2; ++z) // 0 is current block, 1 is nearest along z, y, x
      for(int y = 0; y < 2; ++y)
        for(int x = 0; x < 2; ++x) {
          const Real distx = x? dist_nbix : dist_bix;
          const Real disty = y? dist_nbiy : dist_biy;
          const Real distz = z? dist_nbiz : dist_biz;
          const Real act = action(x? nbix : bix, y? nbiy : biy, z? nbiz : biz);
          const Real dist = std::sqrt(distx*distx + disty*disty + distz*distz);
          const Real weight = std::max( (NB - dist)/NB, (Real) 0);
          weighted_sum_act += act * weight;
          sum_acts_weights += weight;
        }
    //return sum_acts_weights;
    return weighted_sum_act / std::max(sum_acts_weights, (double) 1e-16);
  }

  void set(const double act, const int bix, const int biy, const int biz)
  {
    action(bix, biy, biz) = act;
  }

  const double & action(int bix, int biy, int biz) const
  {
    return actions[(biz+NZ) % NZ][(biy+NY) % NY][(bix+NX) % NX];
  }
  double & action(int bix, int biy, int biz)
  {
    return actions[(biz+NZ) % NZ][(biy+NY) % NY][(bix+NX) % NX];
  }
};

inline Real facFilter(const int i, const int j, const int k)
{
  if (std::abs(i) + std::abs(j) + std::abs(k) == 3)         // Corner cells
    return 1.0/64;
  else if (std::abs(i) + std::abs(j) + std::abs(k) == 2)    // Side-Corner cells
    return 2.0/64;
  else if (std::abs(i) + std::abs(j) + std::abs(k) == 1)    // Side cells
    return 4.0/64;
  else if (std::abs(i) + std::abs(j) + std::abs(k) == 0)    // Center cells
    return 8.0/64;
  else return 0;
}

struct FilteredQuatities
{
  Real u  = 0.0, v  = 0.0, w  = 0.0;
  Real uu = 0.0, uv = 0.0, uw = 0.0;
  Real vv = 0.0, vw = 0.0, ww = 0.0;
  Real S_xx = 0.0, S_xy = 0.0, S_xz = 0.0;
  Real S_yy = 0.0, S_yz = 0.0, S_zz = 0.0;

  FilteredQuatities(Lab& lab, const int ix, const int iy, const int iz, const Real h)
  {
    for (int k=-1; k<2; ++k)
    for (int j=-1; j<2; ++j)
    for (int i=-1; i<2; ++i)
    {
      const Real f = facFilter(i,j,k);
      const auto & L = lab(ix+i, iy+j, iz+k);
      u  += f * L.u;     v  += f * L.v;     w  += f * L.w;
      uu += f * L.u*L.u; uv += f * L.u*L.v; uw += f * L.u*L.w;
      vv += f * L.v*L.v; vw += f * L.v*L.w; ww += f * L.w*L.w;
      const auto & LW=lab(ix+i-1, iy+j, iz+k), & LE=lab(ix+i+1, iy+j, iz+k);
      const auto & LS=lab(ix+i, iy+j-1, iz+k), & LN=lab(ix+i, iy+j+1, iz+k);
      const auto & LF=lab(ix+i, iy+j, iz+k-1), & LB=lab(ix+i, iy+k, iz+k+1);
      S_xx += f * (LE.u-LW.u)             / (2*h);
      S_xy += f * (LE.v-LW.v + LN.u-LS.u) / (4*h);
      S_xz += f * (LE.w-LW.w + LB.u-LF.u) / (4*h);
      S_yy += f * (LN.v-LS.v)             / (2*h);
      S_yz += f * (LN.w-LS.w + LB.v-LF.v) / (4*h);
      S_zz += f * (LB.w-LF.w)             / (2*h);
    }
  }
};

// Product of two symmetric matrices stored as 1D vectors with 6 elts {M_00, M_01, M_02,
//                                                                           M_11, M_12,
//                                                                                 M_22}
// Returns a symmetric matrix.
inline std::array<Real,6> symProd(const std::array<Real,6> & mat1,
                                  const std::array<Real,6> & mat2)
{
  return { mat1[0]*mat2[0] + mat1[1]*mat2[1] + mat1[2]*mat2[2],
           mat1[0]*mat2[1] + mat1[1]*mat2[3] + mat1[2]*mat2[4],
           mat1[0]*mat2[2] + mat1[1]*mat2[4] + mat1[2]*mat2[5],
           mat1[1]*mat2[1] + mat1[3]*mat2[3] + mat1[4]*mat2[4],
           mat1[1]*mat2[2] + mat1[3]*mat2[4] + mat1[4]*mat2[5],
           mat1[2]*mat2[2] + mat1[4]*mat2[4] + mat1[5]*mat2[5] };
}
// Product of two anti symmetric matrices stored as 1D vector with 3 elts (M_01, M_02, M_12)
// Returns a symmetric matrix.
inline std::array<Real,6> antiSymProd(const std::array<Real,3> & mat1,
                                      const std::array<Real,3> & mat2)
{
  return { - mat1[0]*mat2[0] - mat1[1]*mat2[1],
           - mat1[1]*mat2[2],
             mat1[0]*mat2[2],
           - mat1[0]*mat2[0] - mat1[2]*mat2[2],
           - mat1[0]*mat2[1],
           - mat1[1]*mat2[1] - mat1[2]*mat2[2] };
}
// Returns the Tr[mat1*mat2] with mat1 and mat2 symmetric matrices stored as 1D vector.
inline Real traceOfSymProd(const std::array<Real,6> & mat1,
                           const std::array<Real,6> & mat2)
{
  return mat1[0]*mat2[0] +   mat1[3]*mat2[3]  +   mat1[5]*mat2[5]
     + 2*mat1[1]*mat2[1] + 2*mat1[2]*mat2[2]  + 2*mat1[4]*mat2[4];
}

inline std::vector<Real> germanoIdentity(Lab& lab, const Real h,
                                  const int ix, const int iy, const int iz)
{
  FilteredQuatities fq(lab, ix,iy,iz, h);
  const FluidElement & L = lab(ix, iy, iz);
  const Real shear_t = std::sqrt(2*(pow2(fq.S_xx)+pow2(fq.S_yy)+pow2(fq.S_zz))
                               + 4*(pow2(fq.S_xy)+pow2(fq.S_yz)+pow2(fq.S_xz)));
  const Real traceTerm = (fq.uu+fq.vv+fq.ww -fq.u*fq.u -fq.v*fq.v -fq.w*fq.w)/3;
  const Real l_xx = fq.uu - fq.u * fq.u - traceTerm;
  const Real l_xy = fq.uv - fq.u * fq.v;
  const Real l_xz = fq.uw - fq.u * fq.w;
  const Real l_yy = fq.vv - fq.v * fq.v - traceTerm;
  const Real l_yz = fq.vw - fq.v * fq.w;
  const Real l_zz = fq.ww - fq.w * fq.w - traceTerm;
  Real t_xx = - 2 * L.chi * (4*h*h) * shear_t * fq.S_xx;
  Real t_xy = - 2 * L.chi * (4*h*h) * shear_t * fq.S_xy;
  Real t_xz = - 2 * L.chi * (4*h*h) * shear_t * fq.S_xz;
  Real t_yy = - 2 * L.chi * (4*h*h) * shear_t * fq.S_yy;
  Real t_yz = - 2 * L.chi * (4*h*h) * shear_t * fq.S_yz;
  Real t_zz = - 2 * L.chi * (4*h*h) * shear_t * fq.S_zz;
  t_xx -= (t_xx + t_yy + t_zz)/3; // remove trace
  t_yy -= (t_xx + t_yy + t_zz)/3; // remove trace
  t_zz -= (t_xx + t_yy + t_zz)/3; // remove trace

  Real tau_xx = 0, tau_xy = 0, tau_xz = 0, tau_yy = 0, tau_yz = 0, tau_zz = 0;
  for (int i = -1; i < 2; ++i)
  for (int j = -1; j < 2; ++j)
  for (int k = -1; k < 2; ++k)
  {
    const Real f = facFilter(i,j,k);
    //const FluidElement &LL=lab(ix+i,iy+j,iz+k);
    const auto & LW = lab(ix+i-1,iy+j,iz+k), & LE = lab(ix+i+1,iy+j,iz+k);
    const auto & LS = lab(ix+i,iy+j-1,iz+k), & LN = lab(ix+i,iy+j+1,iz+k);
    const auto & LF = lab(ix+i,iy+j,iz+k-1), & LB = lab(ix+i,iy+j,iz+k+1);
    const Real dudx = LE.u-LW.u, dvdx = LE.v-LW.v, dwdx = LE.w-LW.w;
    const Real dudy = LN.u-LS.u, dvdy = LN.v-LS.v, dwdy = LN.w-LS.w;
    const Real dudz = LB.u-LF.u, dvdz = LB.v-LF.v, dwdz = LB.w-LF.w;
    const Real shear = std::sqrt( 2*(dudx*dudx) + 2*(dvdy*dvdy) + 2*(dwdz*dwdz)
                + pow2(dudy+dvdx) + pow2(dudz+dwdx) + pow2(dwdy+dvdz) ) / (2*h);
    tau_xx -= f * 2 * L.chi * shear * h*h *  dudx         / (2*h);
    tau_xy -= f * 2 * L.chi * shear * h*h * (dudy + dvdx) / (4*h);
    tau_xz -= f * 2 * L.chi * shear * h*h * (dudz + dwdx) / (4*h);
    tau_yy -= f * 2 * L.chi * shear * h*h *  dvdy         / (2*h);
    tau_yz -= f * 2 * L.chi * shear * h*h * (dwdy + dvdz) / (4*h);
    tau_zz -= f * 2 * L.chi * shear * h*h *  dwdz         / (2*h);
  }

  return {     ( l_xx - t_xx + tau_xx - (tau_xx + tau_yy + tau_zz)/3),
           2 * ( l_xy - t_xy + tau_xy                               ),
           2 * ( l_xz - t_xz + tau_xz                               ),
               ( l_yy - t_yy + tau_yy - (tau_xx + tau_yy + tau_zz)/3),
           2 * ( l_yz - t_yz + tau_yz                               ),
               ( l_zz - t_zz + tau_zz - (tau_xx + tau_yy + tau_zz)/3) };
}

using rlApi_t = std::function<Real(const std::array<Real, 9> &, const size_t,
                                   const size_t,const int,const int,const int)>;
using locRewF_t = std::function<void(const size_t blockID, Lab & lab)>;

class KernelSGS_RL
{
 private:
  const rlApi_t & sendStateRecvAct;
  const locRewF_t & computeNextLocalRew;
  ActionInterpolator & actInterp;
  const HITstatistics & stats;
  const Real scaleVel, scaleGrad, scaleLap;

  Real sqrtDist(const Real val) const {
    return val>=0? std::sqrt(val) : -std::sqrt(-val);
  };
  Real cbrtDist(const Real val) const {
    return std::cbrt(val);
  };
  Real frthDist(const Real val) const {
    return val>=0? std::sqrt(std::sqrt(val)) : -std::sqrt(std::sqrt(-val));
  };

  std::array<Real,5> popeInvariants(
    const Real d1udx1, const Real d1vdx1, const Real d1wdx1,
    const Real d1udy1, const Real d1vdy1, const Real d1wdy1,
    const Real d1udz1, const Real d1vdz1, const Real d1wdz1) const
  {
    const std::array<Real,6> S = { d1udx1, (d1vdx1 + d1udy1)/2,
      (d1wdx1 + d1udz1)/2, d1vdy1, (d1wdy1 + d1vdz1)/2, d1wdz1 };
    const std::array<Real,3> R = {
      (d1vdx1 - d1udy1)/2, (d1wdx1 - d1udz1)/2, (d1wdy1 - d1vdz1)/2 };
    const std::array<Real,6> S2  = symProd(S, S);
    const std::array<Real,6> R2  = antiSymProd(R, R);
    //const std::vector<Real> R2S = symProd(R2, S);
    return { sqrtDist(S2[0] + S2[3] + S2[5]),   // Tr(S^2)
             sqrtDist(R2[0] + R2[3] + R2[5]),   // Tr(R^2)
             cbrtDist(traceOfSymProd(S2, S)),   // Tr(S^3)
             cbrtDist(traceOfSymProd(R2, S)),   // Tr(R^2.S)
             frthDist(traceOfSymProd(R2,S2)) }; // Tr(R^2.S^2)
  }

  std::array<Real,3> mainMatInvariants(
    const Real xx, const Real xy, const Real xz,
    const Real yx, const Real yy, const Real yz,
    const Real zx, const Real zy, const Real zz) const
  {
    const Real I1 = xx + yy + zz; // Tr(Mat)
    // ( Tr(Mat)^2 - Tr(Mat^2) ) / 2:
    const Real I2 = xx*yy + yy*zz + xx*zz - xy*yx - yz*zy - xz*zx;
    // Det(Mat):
    const Real I3 = xy*yz*zx + xz*yx*zy + xx*yy*zz
                  - xz*yy*zx - xx*yz*zy - xy*yx*zz;
    return {I1, sqrtDist(I2), cbrtDist(I3)};
  }

  template <typename Lab>
  std::array<Real, 9> getState_uniform(Lab& lab, const Real h,
        const int ix, const int iy, const int iz) const
  {
    const Real facGrad = scaleGrad / (2*h), facLap = scaleLap / (h*h);
    const FluidElement &L  = lab(ix, iy, iz);
    const FluidElement &LW = lab(ix - 1, iy, iz), &LE = lab(ix + 1, iy, iz);
    const FluidElement &LS = lab(ix, iy - 1, iz), &LN = lab(ix, iy + 1, iz);
    const FluidElement &LF = lab(ix, iy, iz - 1), &LB = lab(ix, iy, iz + 1);

    const Real d1udx = facGrad*(LE.u-LW.u), d2udx = facLap*(LE.u+LW.u-L.u*2);
    const Real d1vdx = facGrad*(LE.v-LW.v), d2vdx = facLap*(LE.v+LW.v-L.v*2);
    const Real d1wdx = facGrad*(LE.w-LW.w), d2wdx = facLap*(LE.w+LW.w-L.w*2);
    const Real d1udy = facGrad*(LN.u-LS.u), d2udy = facLap*(LN.u+LS.u-L.u*2);
    const Real d1vdy = facGrad*(LN.v-LS.v), d2vdy = facLap*(LN.v+LS.v-L.v*2);
    const Real d1wdy = facGrad*(LN.w-LS.w), d2wdy = facLap*(LN.w+LS.w-L.w*2);
    const Real d1udz = facGrad*(LB.u-LF.u), d2udz = facLap*(LB.u+LF.u-L.u*2);
    const Real d1vdz = facGrad*(LB.v-LF.v), d2vdz = facLap*(LB.v+LF.v-L.v*2);
    const Real d1wdz = facGrad*(LB.w-LF.w), d2wdz = facLap*(LB.w+LF.w-L.w*2);
    const Real S0 = scaleVel * std::sqrt(L.u*L.u + L.v*L.v + L.w*L.w);
    const std::array<double,5> S1 = popeInvariants(d1udx, d1vdx, d1wdx,
                                                   d1udy, d1vdy, d1wdy,
                                                   d1udz, d1vdz, d1wdz);
    const std::array<double,3> S2 = mainMatInvariants(d2udx, d2vdx, d2wdx,
                                                      d2udy, d2vdy, d2wdy,
                                                      d2udz, d2vdz, d2wdz);
    return {S0, S1[0], S1[1], S1[2], S1[3], S1[4], S2[0], S2[1], S2[2]};
  }

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2, 2, 2, false, {FE_U,FE_V,FE_W}};
  //const StencilInfo stencil = StencilInfo(-2,-2,-2, 3,3,3, true, {0,1,2,3});

  KernelSGS_RL(const rlApi_t& api, const locRewF_t& lRew,
        ActionInterpolator& interp, const HITstatistics& _stats,
        const Real _facVel, const Real _facGrad, const Real _facLap) :
        sendStateRecvAct(api), computeNextLocalRew(lRew),
        actInterp(interp), stats(_stats), scaleVel(_facVel),
        scaleGrad(_facGrad), scaleLap(_facLap) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint;
    const size_t thrID = omp_get_thread_num(), blockID = info.blockID;

    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
      const auto state = getState_uniform(lab, h, ix, iy, iz);
      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will have to figure out smth
      // we could compute a local reward here, place as second arg
      o(ix,iy,iz).chi = sendStateRecvAct(state, blockID, thrID, ix,iy,iz);
    }
    //computeNextLocalRew(blockID, lab);
  }

  void state_center(const BlockInfo& info)
  {
    FluidBlock & o = * (FluidBlock *) info.ptrBlock;
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint;
    const size_t thrID = omp_get_thread_num(), blockID = info.blockID;
    const int idx = CUP_BLOCK_SIZE/2 - 1, ipx = CUP_BLOCK_SIZE/2;
    std::array<Real, 9> avgState = {0.};
    const double factor = 1.0 / 8;
    for (int iz = idx; iz <= ipx; ++iz)
    for (int iy = idx; iy <= ipx; ++iy)
    for (int ix = idx; ix <= ipx; ++ix) {
      const auto state = getState_uniform(o, h, ix, iy, iz);
      for (int k = 0; k < 9; ++k) avgState[k] += factor * state[k];
      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will have to figure out smth
      // we could compute a local reward here, place as second arg
    }
    actInterp.set(sendStateRecvAct(avgState, blockID, thrID, idx,idx,idx),
      info.index[0], info.index[1], info.index[2]);
  }

  void apply_actions(const BlockInfo & i) const
  {
    FluidBlock & o = * (FluidBlock *) i.ptrBlock;
    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix)
      o(ix,iy,iz).chi = actInterp(i.index[0], i.index[1], i.index[2], ix,iy,iz);
  }
};

SGS_RL::SGS_RL(SimulationData&s, smarties::Communicator*_comm,
               const int nAgentsPB) : Operator(s), commPtr(_comm),
               nAgentsPerBlock(nAgentsPB)
{
  // TODO relying on chi field does not work is obstacles are present
  // TODO : make sure there are no agents on the same grid point if nAgentsPB>1
  assert(nAgentsPB == 1); // TODO

  std::mt19937& gen = commPtr->getPRNG();
  const std::vector<BlockInfo>& myInfo = sim.vInfo();
  std::uniform_int_distribution<int> disX(0, FluidBlock::sizeX-1);
  std::uniform_int_distribution<int> disY(0, FluidBlock::sizeY-1);
  std::uniform_int_distribution<int> disZ(0, FluidBlock::sizeZ-1);
  agentsIDX.resize(myInfo.size(), -1);
  agentsIDY.resize(myInfo.size(), -1);
  agentsIDZ.resize(myInfo.size(), -1);
  localRewards = std::vector<double>(myInfo.size(), 0);

  for (size_t i=0; i<myInfo.size(); ++i) {
    agentsIDX[i]= disX(gen); agentsIDY[i]= disY(gen); agentsIDZ[i]= disZ(gen);
  }
}

void SGS_RL::run(const double dt, const bool RLinit, const bool RLover,
                 const HITstatistics& stats, const HITtargetData& target,
                 const Real globalR)
{
  sim.startProfiler("SGS_RL");
  smarties::Communicator & comm = * commPtr;
  std::vector<double> nextlocRewards(localRewards.size(), 0);
  ActionInterpolator actInterp( sim.grid->getResidentBlocksPerDimension(2),
                                sim.grid->getResidentBlocksPerDimension(1),
                                sim.grid->getResidentBlocksPerDimension(0) );

  #if 0 // non-dimensionalize wrt flow quantities
    const Real scaleVel = 1 / std::sqrt(stats.tke); // [T/L]
    const Real scaleGrad = stats.tke / stats.dissip_tot; // [T]
    const Real scaleLap = scaleGrad * stats.getKolmogorovL(); // [TL]

    // one element per block is a proper agent: will add seq to train data
    // other are nThreads and are only there for thread safety
    // states get overwritten
    const Real h_nonDim = sim.uniformH() / stats.getKolmogorovL();
    const Real dt_nonDim = dt * stats.dissip_tot / stats.tke;
    const Real tke_nonDim = stats.tke / std::sqrt(stats.dissip_tot * stats.nu);
    const Real visc_nonDim = stats.dissip_visc / stats.dissip_tot;
    const Real lenIn_nonDim = stats.l_integral / stats.lambda;
    const Real inject_nonDim = sim.actualInjectionRate / stats.dissip_tot;
    const Real deltaEn_nonDim = (stats.tke - target.tKinEn) / stats.tke;
  #else // non-dimensionalize wrt *target* flow quantities
    const Real eta = stats.getKolmogorovL(target.epsVis, target.nu);
    const Real scaleVel = 1 / std::sqrt(target.tKinEn); // [T/L]
    const Real scaleGrad = target.tKinEn / target.epsVis; // [T]
    const Real scaleLap = scaleGrad * eta; // [TL]
    const Real h_nonDim = sim.uniformH() / eta;
    const Real dt_nonDim = dt / scaleGrad;
    const Real tke_nonDim = stats.tke / std::sqrt(target.epsVis * target.nu);
    const Real visc_nonDim = stats.dissip_visc / target.epsVis;
    const Real lenIn_nonDim = stats.l_integral / target.lInteg;
    const Real inject_nonDim = stats.dissip_tot / sim.actualInjectionRate;
    const Real deltaEn_nonDim = (stats.tke - target.tKinEn) / target.tKinEn;
  #endif

  std::array<Real, 7> globalS = { h_nonDim, dt_nonDim,
    tke_nonDim, visc_nonDim, lenIn_nonDim, inject_nonDim, deltaEn_nonDim };

  const auto getState = [&] (const std::array<Real,9> & locS) {
    return std::vector<double> { locS[0], locS[1], locS[2], locS[3], locS[4],
      locS[5], locS[6], locS[7], locS[8], globalS[0], globalS[1], globalS[2],
      globalS[3], globalS[4], globalS[5], globalS[6] };
  };

  #if 0 // old setup:
    // Randomly scattered agent-grid-points that sample the policy for Cs.
    // Rest of grid follows the mean of the policy s.t. grad log pi := 0.
    // Therefore only one element per block is a proper agent: will add EP to
    // train data, while other are nThreads agents and are only there for thread
    // safety: their states get overwritten, actions are policy mean.
    const size_t nBlocks = sim.vInfo().size();
    const Uint getAgentID = [&](const size_t blockID, const size_t threadID,
                                const int ix,const int iy,const int iz) {
      const bool bAgent = ix == agentsIDX[blockID] &&
                          iy == agentsIDY[blockID] &&
                          iz == agentsIDZ[blockID];
      return bAgent? blockID : nBlocks + threadID;
    };
  #else // new setup:
    // Agents in block centers and linear interpolate Cs on the grid.
    // The good: (i.) stronger signals for rewards (fewer agents take decisions)
    // (ii.) can use RNN. The bad: Less powerful model, coarse grained state.
    const auto getAgentID = [&](const size_t blockID, const size_t threadID,
                                const int ix,const int iy,const int iz) {
      return blockID;
    };
  #endif

  const rlApi_t Finit = [&](const std::array<Real,9> & locS, const size_t bID,
                   const size_t thrID, const int ix, const int iy, const int iz)
  {
    const size_t agentID = getAgentID(bID, thrID, ix, iy, iz);
    comm.sendInitState(getState(locS), agentID);
    return comm.recvAction(agentID)[0];
  };
  const rlApi_t Fcont = [&](const std::array<Real,9> & locS, const size_t bID,
                   const size_t thrID, const int ix, const int iy, const int iz)
  {
    const size_t agentID = getAgentID(bID, thrID, ix, iy, iz);
    comm.sendState(getState(locS), globalR + localRewards[bID], agentID);
    return comm.recvAction(agentID)[0];
  };
  const rlApi_t Flast = [&](const std::array<Real,9> & locS, const size_t bID,
                   const size_t thrID, const int ix, const int iy, const int iz)
  {
    const size_t agentID = getAgentID(bID, thrID, ix, iy, iz);
    comm.sendLastState(getState(locS), globalR + localRewards[bID], agentID);
    return (Real) 0;
  };
  const rlApi_t sendState = RLinit ? Finit : ( RLover ? Flast : Fcont );

  const locRewF_t computeNextLocalRew = [&] (const size_t bID, Lab& lab)
  {
    const auto ix = agentsIDX[bID], iy = agentsIDY[bID], iz = agentsIDZ[bID];
    const Real h = sim.vInfo()[bID].h_gridpoint;
    const std::vector<Real> germano = germanoIdentity(lab, ix, iy, iz, h);
    nextlocRewards[bID] = -(std::fabs(germano[0])+std::fabs(germano[1]) +
                            std::fabs(germano[2])+std::fabs(germano[3]) +
                            std::fabs(germano[4])+std::fabs(germano[5]))/9;
  };

  KernelSGS_RL K_SGS_RL(sendState, computeNextLocalRew, actInterp,
                        stats, scaleVel, scaleGrad, scaleLap);

  #if 0 // old setup :
    compute<KernelSGS_RL>(K_SGS_RL);
  #else // new setup : (first get actions for block centers, then interpolate)
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < vInfo.size(); ++i) K_SGS_RL.state_center(vInfo[i]);

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < vInfo.size(); ++i) K_SGS_RL.apply_actions(vInfo[i]);
  #endif

  sim.stopProfiler();
  check("SGS_RL");
  localRewards = nextlocRewards;
}

CubismUP_3D_NAMESPACE_END
