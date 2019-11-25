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
    for (int k=-1; k<2; ++k)
    for (int j=-1; j<2; ++j)
    for (int i=-1; i<2; ++i)
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
inline std::array<Real,6> symProd(const std::array<Real,6> & mat1,
                                  const std::array<Real,6> & mat2) {
  std::array<Real,6> ret;
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
inline std::array<Real,6> antiSymProd(const std::array<Real,3> & mat1,
                                      const std::array<Real,3> & mat2)
{
  std::array<Real,6> ret;
  ret[0] = - mat1[0]*mat2[0] - mat1[1]*mat2[1];
  ret[1] = - mat1[1]*mat2[2];
  ret[2] =   mat1[0]*mat2[2];
  ret[3] = - mat1[0]*mat2[0] - mat1[2]*mat2[2];
  ret[4] = - mat1[0]*mat2[1];
  ret[5] = - mat1[1]*mat2[1] - mat1[2]*mat2[2];
  return ret;
}
// Returns the Tr[mat1*mat2] with mat1 and mat2 symmetric matrices stored as 1D vector.
inline Real traceOfSymProd(const std::array<Real,6> & mat1,
                           const std::array<Real,6> & mat2)
{
  Real ret =   mat1[0]*mat2[0] +   mat1[3]*mat2[3]  +   mat1[5]*mat2[5]
           + 2*mat1[1]*mat2[1] + 2*mat1[2]*mat2[2]  + 2*mat1[4]*mat2[4];
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

  ret[0] = (l_xx - (t_xx - tau_xx));
  ret[1] = (l_xy - (t_xy - tau_xy))*2;
  ret[2] = (l_xz - (t_xz - tau_xz))*2;
  ret[3] = (l_yy - (t_yy - tau_yy));
  ret[4] = (l_yz - (t_yz - tau_yz))*2;
  ret[5] = (l_zz - (t_zz - tau_zz));
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
  const Real eps, tke, nu, dt_nonDim;
  const Real invEta = std::pow(eps/(nu*nu*nu), 0.25);
  // const Real scaleL = std::pow(tke, 1.5) / eps; [L]
  const Real scaleVel = 1 / std::sqrt(tke); // [T/L]
  const Real scaleGrad = tke / eps; // [T]
  const Real scaleLap = std::pow(tke, 2.5) / std::pow(eps,2); // [TL]

  Real sqrtDist(const Real val) const {
    return val>=0? std::sqrt(val) : -std::sqrt(-val);
  };
  Real frthDist(const Real val) const {
    return val>=0? std::sqrt(std::sqrt(val)) : -std::sqrt(std::sqrt(-val));
  };

  std::array<Real,5> popeInvariants(
    const Real d1udx1, const Real d1vdx1, const Real d1wdx1,
    const Real d1udy1, const Real d1vdy1, const Real d1wdy1,
    const Real d1udz1, const Real d1vdz1, const Real d1wdz1) const
  {
    const std::array<Real,6> S = {
      d1udx1, (d1vdx1 + d1udy1)/2, (d1wdx1 + d1udz1)/2,
      d1vdy1, (d1wdy1 + d1vdz1)/2, d1wdz1 };

    const std::array<Real,3> R = {
      (d1vdx1 - d1udy1)/2, (d1wdx1 - d1udz1)/2, (d1wdy1 - d1vdz1)/2};

    const std::array<Real,6> S2  = symProd(S, S);
    const std::array<Real,6> R2  = antiSymProd(R, R);
    //const std::vector<Real> R2S = symProd(R2, S);
    std::array<Real,5> ret;

    ret[0] = sqrtDist(S2[0] + S2[3] + S2[5]);  // Tr(S^2)
    ret[1] = sqrtDist(R2[0] + R2[3] + R2[5]);  // Tr(R^2)
    ret[2] = std::cbrt(traceOfSymProd(S2, S)); // Tr(S^3)
    ret[3] = std::cbrt(traceOfSymProd(R2, S)); // Tr(R^2.S)
    ret[4] = frthDist(traceOfSymProd(R2, S2)); // Tr(R^2.S^2)

    return ret;
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
    return {I1, sqrtDist(I2), std::cbrt(I3)};
  }

  std::vector<double> getState_uniform(Lab& lab, const Real h,
        const Real h_nonDim, const int ix, const int iy, const int iz) const
  {
    const Real facGrad = scaleGrad / (2*h), facLap = scaleLap / (h*h);
    const FluidElement &L  = lab(ix, iy, iz);
    const FluidElement &LW = lab(ix - 1, iy, iz), &LE = lab(ix + 1, iy, iz);
    const FluidElement &LS = lab(ix, iy - 1, iz), &LN = lab(ix, iy + 1, iz);
    const FluidElement &LF = lab(ix, iy, iz - 1), &LB = lab(ix, iy, iz + 1);

    const Real d1udx = facGrad*(LE.u-LW.u), d2udx = facLap*(LN.u+LS.u-L.u*6);
    const Real d1vdx = facGrad*(LE.v-LW.v), d2vdx = facLap*(LN.v+LS.v-L.v*6);
    const Real d1wdx = facGrad*(LE.w-LW.w), d2wdx = facLap*(LN.w+LS.w-L.w*6);
    const Real d1udy = facGrad*(LN.u-LS.u), d2udy = facLap*(LE.u+LW.u-L.u*2);
    const Real d1vdy = facGrad*(LN.v-LS.v), d2vdy = facLap*(LE.v+LW.v-L.v*2);
    const Real d1wdy = facGrad*(LN.w-LS.w), d2wdy = facLap*(LE.w+LW.w-L.w*2);
    const Real d1udz = facGrad*(LB.u-LF.u), d2udz = facLap*(LF.u+LB.u-L.u*2);
    const Real d1vdz = facGrad*(LB.v-LF.v), d2vdz = facLap*(LF.v+LB.v-L.v*2);
    const Real d1wdz = facGrad*(LB.w-LF.w), d2wdz = facLap*(LF.w+LB.w-L.w*2);
    const Real S0 = scaleVel * std::sqrt(L.u*L.u + L.v*L.v + L.w*L.w);
    const std::array<double,5> S1 = popeInvariants(d1udx, d1vdx, d1wdx,
                                                  d1udy, d1vdy, d1wdy,
                                                  d1udz, d1vdz, d1wdz);
    const std::array<double,3> S2 = mainMatInvariants(d2udx, d2vdx, d2wdx,
                                                     d2udy, d2vdy, d2wdy,
                                                     d2udz, d2vdz, d2wdz);
    return {dt_nonDim, h_nonDim, S0,
            S1[0], S1[1], S1[2], S1[3], S1[4], S2[0], S2[1], S2[2]};
  }

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2, 2, 2, false, {FE_U,FE_V,FE_W}};
  //const StencilInfo stencil = StencilInfo(-2,-2,-2, 3,3,3, true, {0,1,2,3});

  KernelSGS_RL(const rlApi_t& api, const locRewF_t& lRew, Real _eps, Real _tke,
        Real _nu, Real _dt) : sendStateRecvAct(api), computeNextLocalRew(lRew),
        eps(_eps), tke(_tke), nu(_nu), dt_nonDim(_dt * eps / tke) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint, h_nonDim = h * invEta;
    const size_t thrID = omp_get_thread_num(), blockID = info.blockID;

    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
      const auto state = getState_uniform(lab, h, h_nonDim, ix, iy, iz);
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
                 const Real eps, const Real tke, const Real collectiveReward)
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

  const KernelSGS_RL K_SGS_RL(sendState,computeNextLocalRew,eps,tke,sim.nu,dt);

  compute<KernelSGS_RL>(K_SGS_RL);
  sim.stopProfiler();
  check("SGS_RL");
  localRewards = nextlocRewards;
}

CubismUP_3D_NAMESPACE_END
