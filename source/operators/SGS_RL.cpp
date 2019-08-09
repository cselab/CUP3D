//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/Operator.h"
#include "operators/SGS_RL.h"
#include "Communicators/Communicator.h"

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

inline int getAgentId(const int idx, const int idy, const int idz,
               const std::vector<int> trackedAgentsX,
               const std::vector<int> trackedAgentsY,
               const std::vector<int> trackedAgentsZ)
{
  const int nAgentsPerBlock = trackedAgentsX.size();
  for (int i=0; i<nAgentsPerBlock; ++i)
  {
    if (idx==trackedAgentsX[i] and idy==trackedAgentsY[i] and idz==trackedAgentsZ[i])
      return i;
  }
  return -1;
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

class KernelSGS_RL
{
 private:
  smarties::Communicator& comm;
  const int step;
  const bool timeOut;
  const double reward;
  const size_t nBlocks;
  const size_t nAgentsPerBlock;
  const Real scaleGrads;

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2, 2, 2, false, {FE_U,FE_V,FE_W}};

  KernelSGS_RL(smarties::Communicator& _comm, const int _step,
               const bool _timeOut, const double _rew, const Real _scaleGrads,
               const size_t _nBlocks, const size_t _nAgentsPerBlock) :
  comm(_comm), step(_step), timeOut(_timeOut), reward(_rew), nBlocks(_nBlocks),
  nAgentsPerBlock(_nAgentsPerBlock), scaleGrads(_scaleGrads) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint;
    const int thrID = omp_get_thread_num();
    const size_t nAgents = nAgentsPerBlock * nBlocks;
    using send_t = std::function<void(const std::vector<double>,double,size_t)>;
    const send_t Finit = [&](const std::vector<double> S, double R, size_t ID) {
            comm.sendInitState(S, ID); };
    const send_t Fcont = [&](const std::vector<double> S, double R, size_t ID) {
            comm.sendState(S, R, ID); };
    const send_t Flast = [&](const std::vector<double> S, double R, size_t ID) {
            comm.sendLastState(S, R, ID); };
    const send_t sendState = step == 0 ? Finit : ( timeOut ? Flast : Fcont );

    // Deal with the real agents first
    for (size_t k = 0; k < nAgentsPerBlock; ++k)
    {
      const int ix = o.iAgentX[k], iy = o.iAgentY[k], iz = o.iAgentZ[k];
      const size_t agentID = nAgentsPerBlock*info.blockID + k;
      const auto state = getState_uniform(lab, ix, iy, iz, h, scaleGrads);
      sendState(state, reward, agentID);
      if (!timeOut) o(ix,iy,iz).chi = comm.recvAction(agentID)[0];
    }

    // Then the fake agents
    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix)
    {
      // one element per block is a proper agent: will add seq to train data
      // other are nThreads and are only there for thread safety
      // states get overwritten

      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will have to figure out smth
      const int locAgentID= getAgentId(ix,iy,iz, o.iAgentX,o.iAgentY,o.iAgentZ);

      //std::cout<<"Local Agent id"<< localAgentID << std::endl;
      if (locAgentID>=0) continue;
      const size_t agentID =  nAgents + thrID;
      const auto state = getState_uniform(lab, ix, iy, iz, h, scaleGrads);
      sendState(state, reward, agentID);
      if (!timeOut) o(ix,iy,iz).chi = comm.recvAction(agentID)[0];
    }
  }
};

SGS_RL::SGS_RL(SimulationData& s, smarties::Communicator* _comm,
               const int _step, const bool _timeOut,
               const double _reward, const double _scaleGrads,
               const int _nAgentsPerBlock) : Operator(s), comm(_comm),
               step(_step), timeOut(_timeOut), reward(_reward),
               scaleGrads(_scaleGrads), nAgentsPerBlock(_nAgentsPerBlock)
{}

void SGS_RL::operator()(const double dt)
{
  sim.startProfiler("SGS_RL");
  const KernelSGS_RL K_SGS_RL(*comm, step, timeOut, reward, scaleGrads,
                              sim.vInfo().size(), nAgentsPerBlock);

  compute<KernelSGS_RL>(K_SGS_RL);
  sim.stopProfiler();
  check("SGS_RL");
}

CubismUP_3D_NAMESPACE_END
