//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/SGS_RL.h"
#include "Communicator.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class KernelSGS_RL_nonUniform
{
  private:
  const Real time, mu;
  const Real* const uInf;
  Communicator& comm;
  const bool timeOut;
  const double reward;
  const int nBlocks;
  const envInfo rlSeqInfo = time<=0? INIT_COMM:(timeOut? TRNC_COMM : CONT_COMM);

  public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);

  KernelSGS_RL_nonUniform(double t, double m, Real*const u, Communicator& c,
    bool over, double rew, int nblocks) : time(t), mu(m), uInf(u), comm(c),
    timeOut(over), reward(rew), nBlocks(nblocks) { }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    //const BlkCoeffX & c1x = o.fd_cx.first, & c2x = o.fd_cx.second;
    //const BlkCoeffY & c1y = o.fd_cy.first, & c2y = o.fd_cy.second;
    //const BlkCoeffZ & c1z = o.fd_cz.first, & c2z = o.fd_cz.second;
    const int thrID = omp_get_thread_num();
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // one element per block is a proper agent: will add seq to train data
      // other are nThreads and are only there for thread safety
      // states get overwritten: TODO check that it works properly
      const int agentID = iz==0 && iy==0 && iz==0? info.blockID : nBlocks+thrID;
      const std::vector<double> dummy(27, 0);
      comm.sendState(agentID, rlSeqInfo, dummy, reward);
      /*
        const FluidElement &L =lab(ix,iy,iz);
        const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
        const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
        const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
        const Real d1udx1 = __FD_2ND(ix, c1x, LW.u, L.u, LE.u);
        const Real d1udy1 = __FD_2ND(iy, c1y, LS.u, L.u, LN.u);
        const Real d1udz1 = __FD_2ND(iz, c1z, LF.u, L.u, LB.u);
        const Real d1vdx1 = __FD_2ND(ix, c1x, LW.v, L.v, LE.v);
        const Real d1vdy1 = __FD_2ND(iy, c1y, LS.v, L.v, LN.v);
        const Real d1vdz1 = __FD_2ND(iz, c1z, LF.v, L.v, LB.v);
        const Real d1wdx1 = __FD_2ND(ix, c1x, LW.w, L.w, LE.w);
        const Real d1wdy1 = __FD_2ND(iy, c1y, LS.w, L.w, LN.w);
        const Real d1wdz1 = __FD_2ND(iz, c1z, LF.w, L.w, LB.w);
        const Real d2udx2 = __FD_2ND(ix, c2x, LW.u, L.u, LE.u);
        const Real d2udy2 = __FD_2ND(iy, c2y, LS.u, L.u, LN.u);
        const Real d2udz2 = __FD_2ND(iz, c2z, LF.u, L.u, LB.u);
        const Real d2vdx2 = __FD_2ND(ix, c2x, LW.v, L.v, LE.v);
        const Real d2vdy2 = __FD_2ND(iy, c2y, LS.v, L.v, LN.v);
        const Real d2vdz2 = __FD_2ND(iz, c2z, LF.v, L.v, LB.v);
        const Real d2wdx2 = __FD_2ND(ix, c2x, LW.w, L.w, LE.w);
        const Real d2wdy2 = __FD_2ND(iy, c2y, LS.w, L.w, LN.w);
        const Real d2wdz2 = __FD_2ND(iz, c2z, LF.w, L.w, LB.w);
        const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
        const Real duD = d2udx2 + d2udy2 + d2udz2;
        const Real dvD = d2vdx2 + d2vdy2 + d2vdz2;
        const Real dwD = d2wdx2 + d2wdy2 + d2wdz2;
        const Real duA = u * d1udx1 + v * d1udy1 + w * d1udz1;
        const Real dvA = u * d1vdx1 + v * d1vdy1 + w * d1vdz1;
        const Real dwA = u * d1wdx1 + v * d1wdy1 + w * d1wdz1;
      */
      std::vector<double> act = comm.recvAction(agentID);
      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will ahve to figure out smth
      o(ix,iy,iz).chi = act[0]; /* some function */
    }
  }
};

SGS_RL::SGS_RL(SimulationData & s, Communicator* c, bool _timeOut, double rew) :
  Operator(s), comm(c), timeOut(_timeOut), reward(rew) { }

void SGS_RL::operator()(const double dt)
{
  sim.startProfiler("SGS_RL Kernel");
  const KernelSGS_RL_nonUniform K(sim.time, sim.nu, sim.uinf.data(), * comm,
    timeOut, reward, sim.vInfo().size());
  compute<KernelSGS_RL_nonUniform>(K);
  sim.stopProfiler();
  check("SGS_RL");
}

CubismUP_3D_NAMESPACE_END
