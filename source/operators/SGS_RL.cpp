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

#include <functional>
CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

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

struct filteredQuatities
{
  Real u  = 0., v  = 0., w  = 0.;
  Real uu = 0., uv = 0., uw = 0.;
  Real vv = 0., vw = 0., ww = 0.;

  Real S_xx = 0., S_xy = 0., S_xz = 0.;
  Real S_yy = 0., S_yz = 0., S_zz = 0.;

  filteredQuatities(Lab& lab, const int ix, const int iy, const int iz, const Real h)
  {
    for (int i=-1; i<2; i++)
    for (int j=-1; j<2; j++)
    for (int k=-1; k<2; k++)
    {
      const Real f = facFilter(i,j,k);
      const FluidElement &L =lab(ix+i, iy+j, iz+k);
      u  += f * L.u;     v  += f * L.v;     w  += f * L.w;
      uu += f * L.u*L.u; uv += f * L.u*L.v; uw += f * L.u*L.w;
      vv += f * L.v*L.v; vw += f * L.v*L.w; ww += f * L.w*L.w;

      const FluidElement &LW=lab(ix+i-1, iy+j,   iz+k  );
      const FluidElement &LE=lab(ix+i+1, iy+j,   iz+k  );
      const FluidElement &LS=lab(ix+i,   iy+j-1, iz+k  );
      const FluidElement &LN=lab(ix+i,   iy+j+1, iz+k  );
      const FluidElement &LF=lab(ix+i,   iy+j,   iz+k-1);
      const FluidElement &LB=lab(ix+i,   iy+k,   iz+k+1);

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
inline std::vector<Real> symProd(const std::vector<Real> mat1,
                                 const std::vector<Real> mat2)
{
  assert(mat1.size()==6 && mat2.size()==6);
  std::vector<Real> ret(6, 0);
  ret[0] = mat1[0]*mat2[0] + mat1[1]*mat2[1] + mat1[2]*mat2[2];
  ret[1] = mat1[0]*mat2[1] + mat1[1]*mat2[3] + mat1[2]*mat2[4];
  ret[2] = mat1[0]*mat2[2] + mat1[1]*mat2[4] + mat1[2]*mat2[5];
  ret[3] = mat1[1]*mat2[1] + mat1[3]*mat2[3] + mat1[4]*mat2[4];
  ret[4] = mat1[1]*mat2[2] + mat1[3]*mat2[4] + mat1[4]*mat2[5];
  ret[5] = mat1[2]*mat2[2] + mat1[4]*mat2[4] + mat1[5]*mat2[5];
  return ret;
}
// Product of two anti symmetric matrices stored as 1D vector with 3 elts (M_01, M_02, M_12)
// Returns a symmetric matrix.
inline std::vector<Real> antiSymProd(const std::vector<Real> mat1,
                                     const std::vector<Real> mat2)
{
  assert(mat1.size()==3 && mat2.size()==3);
  std::vector<Real> ret(6, 0);
  ret[0] = - mat1[0]*mat2[0] - mat1[1]*mat2[1];
  ret[1] = - mat1[1]*mat2[2];
  ret[2] =   mat1[0]*mat2[2];
  ret[3] = - mat1[0]*mat2[0] - mat1[2]*mat2[2];
  ret[4] = - mat1[0]*mat2[1];
  ret[5] = - mat1[1]*mat2[1] - mat1[2]*mat2[2];
  return ret;
}
// Returns the Tr[mat1*mat2] with mat1 and mat2 symmetric matrices stored as 1D vector.
inline Real traceOfProd(const std::vector<Real> mat1,
                        const std::vector<Real> mat2)
{
  assert(mat1.size()==6 && mat2.size()==6);
  Real ret =   mat1[0]*mat2[0] +   mat1[3]*mat2[3]  +   mat1[5]*mat2[5]
           + 2*mat1[1]*mat2[1] + 2*mat1[2]*mat2[2]  + 2*mat1[4]*mat2[4];
  return ret;
}

inline std::vector<Real> flowInvariants(
  const Real d1udx1, const Real d1vdx1, const Real d1wdx1,
  const Real d1udy1, const Real d1vdy1, const Real d1wdy1,
  const Real d1udz1, const Real d1vdz1, const Real d1wdz1)
{
  const std::vector<Real> S = {
    d1udx1, (d1vdx1 + d1udy1)/2, (d1wdx1 + d1udz1)/2,
    d1vdy1, (d1wdy1 + d1vdz1)/2, d1wdz1 };

  const std::vector<Real> R = {
    (d1vdx1 - d1udy1)/2, (d1wdx1 - d1udz1)/2, (d1wdy1 - d1vdz1)/2};

  const std::vector<Real> S2  = symProd(S, S);
  const std::vector<Real> R2  = antiSymProd(R, R);
  const std::vector<Real> R2S = symProd(R2, S);
  std::vector<Real> ret(5, 0);
  ret[0] = S2[0] + S2[3] + S2[5]; // Tr(S^2)
  ret[1] = R2[0] + R2[3] + R2[5]; // Tr(R^2)
  ret[2] = traceOfProd(S2, S);    // Tr(S^3)
  ret[3] = traceOfProd(R2, S);    // Tr(R^2.S)
  ret[4] = traceOfProd(R2, S2);   // Tr(R^2.S^2)
  return ret;
}

inline std::vector<double> getState_uniform(Lab& lab,
                                     const int ix, const int iy, const int iz,
                                     const Real h, const Real scaleGrads)
{
  //const FluidElement &L  = lab(ix, iy, iz);
  const FluidElement &LW = lab(ix - 1, iy, iz),
                     &LE = lab(ix + 1, iy, iz);
  const FluidElement &LS = lab(ix, iy - 1, iz),
                     &LN = lab(ix, iy + 1, iz);
  const FluidElement &LF = lab(ix, iy, iz - 1),
                     &LB = lab(ix, iy, iz + 1);

  const Real d1udx1= LE.u-LW.u, d1vdx1= LE.v-LW.v, d1wdx1= LE.w-LW.w;
  const Real d1udy1= LN.u-LS.u, d1vdy1= LN.v-LS.v, d1wdy1= LN.w-LS.w;
  const Real d1udz1= LB.u-LF.u, d1vdz1= LB.v-LF.v, d1wdz1= LB.w-LF.w;
  const Real fac = scaleGrads / (2*h);
  #ifdef SGSRL_STATE_INVARIANTS
  const std::vector<double> ret =
    flowInvariants(d1udx1 * fac, d1vdx1 * fac, d1wdx1 * fac,
                   d1udy1 * fac, d1vdy1 * fac, d1wdy1 * fac,
                   d1udz1 * fac, d1vdz1 * fac, d1wdz1 * fac);
  #else
  const std::vector<double> ret = {d1udx1 * fac, d1vdx1 * fac, d1wdx1 * fac,
                                   d1udy1 * fac, d1vdy1 * fac, d1wdy1 * fac,
                                   d1udz1 * fac, d1vdz1 * fac, d1wdz1 * fac};
  #endif
  return ret;
}



inline std::vector<Real> germanoIdentity(Lab& lab, const Real h,
                                  const int ix, const int iy, const int iz)
{
  std::vector<Real> ret(6,0);

  filteredQuatities fq(lab,ix,iy,iz, h);
  const FluidElement &L =lab(ix, iy, iz);

  const Real shear_t = std::sqrt(
      2 * (fq.S_xx * fq.S_xx + fq.S_yy * fq.S_yy + fq.S_zz * fq.S_zz) +
      4 * (fq.S_xy * fq.S_xy + fq.S_yz * fq.S_yz + fq.S_xz * fq.S_xz));

  Real traceTerm = (fq.uu+fq.vv+fq.ww -fq.u*fq.u -fq.v*fq.v -fq.w*fq.w)/3;
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
  traceTerm = (t_xx + t_yy + t_zz)/3;
  t_xx -= traceTerm;
  t_yy -= traceTerm;
  t_zz -= traceTerm;


  Real tau_xx = 0., tau_xy = 0., tau_xz = 0.;
  Real tau_yy = 0., tau_yz = 0., tau_zz = 0.;

  for (int i=-1; i<2; i++)
  for (int j=-1; j<2; j++)
  for (int k=-1; k<2; k++)
  {
    const Real f = facFilter(i,j,k);
    //const FluidElement &LL=lab(ix+i,iy+j,iz+k);
    const FluidElement &LW=lab(ix+i-1,iy+j,iz+k), &LE=lab(ix+i+1,iy+j,iz+k);
    const FluidElement &LS=lab(ix+i,iy+j-1,iz+k), &LN=lab(ix+i,iy+j+1,iz+k);
    const FluidElement &LF=lab(ix+i,iy+j,iz+k-1), &LB=lab(ix+i,iy+j,iz+k+1);
    const Real dudx = LE.u-LW.u, dvdx = LE.v-LW.v, dwdx = LE.w-LW.w;
    const Real dudy = LN.u-LS.u, dvdy = LN.v-LS.v, dwdy = LN.w-LS.w;
    const Real dudz = LB.u-LF.u, dvdz = LB.v-LF.v, dwdz = LB.w-LF.w;

    const Real shear = std::sqrt(
       2*(dudx*dudx) + 2*(dvdy*dvdy) + 2*(dwdz*dwdz)
       + (dudy+dvdx)*(dudy+dvdx)
       + (dudz+dwdx)*(dudz+dwdx)
       + (dwdy+dvdz)*(dwdy+dvdz) ) / (2*h);

    tau_xx -= f * 2 * L.chi * shear * h*h *  dudx         / (2*h);
    tau_xy -= f * 2 * L.chi * shear * h*h * (dudy + dvdx) / (4*h);
    tau_xz -= f * 2 * L.chi * shear * h*h * (dudz + dwdx) / (4*h);
    tau_yy -= f * 2 * L.chi * shear * h*h *  dvdy         / (2*h);
    tau_yz -= f * 2 * L.chi * shear * h*h * (dwdy + dvdz) / (4*h);
    tau_zz -= f * 2 * L.chi * shear * h*h *  dwdz         / (2*h);
  }
  traceTerm = (tau_xx + tau_yy + tau_zz)/3;
  tau_xx -= traceTerm;
  tau_yy -= traceTerm;
  tau_zz -= traceTerm;

  ret[0] = (l_xx - (t_xx - tau_xx))/l_xx;
  ret[1] = (l_xy - (t_xy - tau_xy))/l_xy*2;
  ret[2] = (l_xz - (t_xz - tau_xz))/l_xz*2;
  ret[3] = (l_yy - (t_yy - tau_yy))/l_yy;
  ret[4] = (l_yz - (t_yz - tau_yz))/l_yz*2;
  ret[5] = (l_zz - (t_zz - tau_zz))/l_zz;
  return ret;
}

using rlApi_t = std::function<Real(const std::vector<double>, const double,
                                   const size_t, const size_t,
                                   const int, const int, const int)>;
using locRewF_t = std::function<void(const size_t blockID, Lab & lab)>;

class KernelSGS_RL
{
 private:
  const rlApi_t& sendStateRecvAct;
  const locRewF_t& computeNextLocalRew;
  const Real scaleGrads;

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  //const StencilInfo stencil{-1,-1,-1, 2, 2, 2, false, {FE_U,FE_V,FE_W}};
  const StencilInfo stencil = StencilInfo(-2,-2,-2, 3,3,3, true, {0,1,2,3});

  KernelSGS_RL(const rlApi_t& api, const locRewF_t& lRew, const Real _scaleG) :
    sendStateRecvAct(api), computeNextLocalRew(lRew), scaleGrads(_scaleG) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint;
    const size_t thrID = omp_get_thread_num(), blockID = info.blockID;

    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
      const auto state = getState_uniform(lab, ix, iy, iz, h, scaleGrads);
      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will have to figure out smth
      // we could compute a local reward here, place as second arg
      o(ix,iy,iz).chi = sendStateRecvAct(state, 0, blockID, thrID, ix,iy,iz);
    }
    //computeNextLocalRew(blockID, lab);
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
  std::uniform_int_distribution<int> distX(0, FluidBlock::sizeX-1);
  std::uniform_int_distribution<int> distY(0, FluidBlock::sizeY-1);
  std::uniform_int_distribution<int> distZ(0, FluidBlock::sizeZ-1);
  agentsIDX.resize(myInfo.size(), -1);
  agentsIDY.resize(myInfo.size(), -1);
  agentsIDZ.resize(myInfo.size(), -1);
  localRewards = std::vector<double>(myInfo.size(), 0);

  for (size_t i=0; i<myInfo.size(); ++i) {
    agentsIDX[i] = distX(gen);
    agentsIDY[i] = distY(gen);
    agentsIDZ[i] = distZ(gen);
  }
}

void SGS_RL::run(const double dt, const bool RLinit, const bool RLover,
                 const Real stateScaling, const Real collectiveReward)
{
  sim.startProfiler("SGS_RL");
  smarties::Communicator & comm = * commPtr;
  const size_t nBlocks = sim.vInfo().size();
  std::vector<double> nextlocRewards(localRewards.size(), 0);
  // one element per block is a proper agent: will add seq to train data
  // other are nThreads and are only there for thread safety
  // states get overwritten

  const rlApi_t Finit = [&](const std::vector<double>&S, const double localRew,
                            const size_t blockID, const size_t threadID,
                            const int ix,const int iy,const int iz)
  {
    const bool bAgent = ix == agentsIDX[blockID] &&
                        iy == agentsIDY[blockID] &&
                        iz == agentsIDZ[blockID];
    const size_t agentID = bAgent? blockID : nBlocks + threadID;
    comm.sendInitState(S, agentID);
    return comm.recvAction(agentID)[0];
  };
  const rlApi_t Fcont = [&](const std::vector<double>&S, const double localRew,
                            const size_t blockID, const size_t threadID,
                            const int ix,const int iy,const int iz)
  {
    const bool bAgent = ix == agentsIDX[blockID] &&
                        iy == agentsIDY[blockID] &&
                        iz == agentsIDZ[blockID];
    const size_t agentID = bAgent? blockID : nBlocks + threadID;
    const Real R = collectiveReward + localRewards[blockID];
    comm.sendState(S, R, agentID);
    return comm.recvAction(agentID)[0];
  };
  const rlApi_t Flast = [&](const std::vector<double>&S, const double localRew,
                            const size_t blockID, const size_t threadID,
                            const int ix,const int iy,const int iz)
  {
    const bool bAgent = ix == agentsIDX[blockID] &&
                        iy == agentsIDY[blockID] &&
                        iz == agentsIDZ[blockID];
    const size_t agentID = bAgent? blockID : nBlocks + threadID;
    const Real R = collectiveReward + localRewards[blockID];
    comm.sendLastState(S, R, agentID);
    return (Real) 0;
  };
  const rlApi_t sendState = RLinit ? Finit : ( RLover ? Flast : Fcont );

  const locRewF_t computeNextLocalRew = [&] (const size_t blockID, Lab& lab)
  {
    const auto ix = agentsIDX[blockID];
    const auto iy = agentsIDY[blockID];
    const auto iz = agentsIDZ[blockID];
    const Real h = sim.vInfo()[blockID].h_gridpoint;
    const std::vector<Real> germano = germanoIdentity(lab, ix, iy, iz, h);
    nextlocRewards[blockID] = -(std::fabs(germano[0])+std::fabs(germano[1]) +
                                std::fabs(germano[2])+std::fabs(germano[3]) +
                                std::fabs(germano[4])+std::fabs(germano[5]))/9;
  };

  const KernelSGS_RL K_SGS_RL(sendState, computeNextLocalRew, stateScaling);

  compute<KernelSGS_RL>(K_SGS_RL);
  sim.stopProfiler();
  check("SGS_RL");
  localRewards = nextlocRewards;
}

CubismUP_3D_NAMESPACE_END
