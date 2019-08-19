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

using rlApi_t = std::function<Real(const std::vector<double>, const double,
                                   const size_t, const size_t,
                                   const int, const int, const int)>;
class KernelSGS_RL
{
 private:
  const rlApi_t& sendStateRecvAct;
  const Real scaleGrads;

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2, 2, 2, false, {FE_U,FE_V,FE_W}};

  KernelSGS_RL(const rlApi_t& api, const Real _scaleGrads) :
    sendStateRecvAct(api), scaleGrads(_scaleGrads) {}

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
    const Real R = collectiveReward; // can weigh with local
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
    const Real R = collectiveReward; // can weigh with local
    comm.sendLastState(S, R, agentID);
    return (Real) 0;
  };
  const rlApi_t sendState = RLinit ? Finit : ( RLover ? Flast : Fcont );

  const KernelSGS_RL K_SGS_RL(sendState, stateScaling);

  compute<KernelSGS_RL>(K_SGS_RL);
  sim.stopProfiler();
  check("SGS_RL");
}

CubismUP_3D_NAMESPACE_END
