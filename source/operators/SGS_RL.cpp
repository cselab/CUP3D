//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/Operator.h"
#include "operators/SGS_RL.h"
#include "Communicator.h"

CubismUP_3D_NAMESPACE_BEGIN using namespace cubism;



// Product of two symmetric matrices stored as 1D vectors with 6 elts {M_00, M_01, M_02,
//                                                                           M_11, M_12,
//                                                                                 M_22}
// Returns a symmetric matrix.
std::vector<Real> symProd(const std::vector<Real> mat1, const std::vector<Real> mat2)
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
std::vector<Real> antiSymProd(const std::vector<Real> mat1, const std::vector<Real> mat2)
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
Real traceOfProd(const std::vector<Real> mat1, const std::vector<Real> mat2)
{
  assert(mat1.size()==6 && mat2.size()==6);
  Real ret =   mat1[0]*mat2[0] +   mat1[3]*mat2[3]  +   mat1[5]*mat2[5]
           + 2*mat1[1]*mat2[1] + 2*mat1[2]*mat2[2]  + 2*mat1[4]*mat2[4];
  return ret;
}

std::vector<Real> flowInvariants(const Real d1udx1, const Real d1vdx1, const Real d1wdx1,
                                 const Real d1udy1, const Real d1vdy1, const Real d1wdy1,
                                 const Real d1udz1, const Real d1vdz1, const Real d1wdz1)
{
  const std::vector<Real> S = {d1udx1, 0.5*(d1vdx1 + d1udy1), 0.5*(d1wdx1 + d1udz1),
                                                d1vdy1      , 0.5*(d1wdy1 + d1vdz1),
                                                                      d1wdz1      };

  const std::vector<Real> R = {0.5*(d1vdx1 - d1udy1), 0.5*(d1wdx1 - d1udz1), 0.5*(d1wdy1 - d1vdz1)};

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

int getAgentId(const int idx, const int idy, const int idz,
               const std::vector<int> trackedAgentsX,
               const std::vector<int> trackedAgentsY,
               const std::vector<int> trackedAgentsZ){
  const int nAgentsPerBlock = trackedAgentsX.size();
  for (int i=0; i<nAgentsPerBlock; i++)
  {
    if (idx==trackedAgentsX[i] and idy==trackedAgentsY[i] and idz==trackedAgentsZ[i])
      return i;
  }
  return -1;
}

std::vector<double> getState_uniform(Lab& lab, const int ix, const int iy, const int iz, const Real h){
  const FluidElement &L  = lab(ix, iy, iz);
  const FluidElement &LW = lab(ix - 1, iy, iz),
                     &LE = lab(ix + 1, iy, iz);
  const FluidElement &LS = lab(ix, iy - 1, iz),
                     &LN = lab(ix, iy + 1, iz);
  const FluidElement &LF = lab(ix, iy, iz - 1),
                     &LB = lab(ix, iy, iz + 1);

  const Real d1udx1= LE.u-LW.u, d1vdx1= LE.v-LW.v, d1wdx1= LE.w-LW.w;
  const Real d1udy1= LN.u-LS.u, d1vdy1= LN.v-LS.v, d1wdy1= LN.w-LS.w;
  const Real d1udz1= LB.u-LF.u, d1vdz1= LB.v-LF.v, d1wdz1= LB.w-LF.w;
  /*
  const std::vector<double> ret = flowInvariants(d1udx1/(2*h), d1vdx1/(2*h), d1wdx1/(2*h),
                                                 d1udy1/(2*h), d1vdy1/(2*h), d1wdy1/(2*h),
                                                 d1udz1/(2*h), d1vdz1/(2*h), d1wdz1/(2*h));
  */
  const std::vector<double> ret = {d1udx1/(2*h), d1vdx1/(2*h), d1wdx1/(2*h),
                                   d1udy1/(2*h), d1vdy1/(2*h), d1wdy1/(2*h),
                                   d1udz1/(2*h), d1vdz1/(2*h), d1wdz1/(2*h)};

  return ret;
}

class KernelSGS_RL {
 private:
  Communicator& comm;
  const int step;
  const bool timeOut;
  const bool evalStep;
  const double reward;
  const size_t nBlocks;
  const int nAgentsPerBlock;
  const envInfo rlSeqInfo =
      (step == 0) ? INIT_COMM : (timeOut ? TRNC_COMM : CONT_COMM);

 public:
  const std::array<int, 3> stencil_start = {-1, -1, -1},
                           stencil_end = {2, 2, 2};
  const StencilInfo stencil =
      StencilInfo(-1, -1, -1, 2, 2, 2, false, 3, 1, 2, 3);

  KernelSGS_RL(Communicator& _comm, const int _step, const bool _timeOut,
               const bool _evalStep, const double _reward, const size_t _nBlocks,
               const int _nAgentsPerBlock)
      : comm(_comm),
        step(_step),
        timeOut(_timeOut),
        evalStep(_evalStep),
        reward(_reward),
        nBlocks(_nBlocks),
        nAgentsPerBlock(_nAgentsPerBlock) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab& lab, const BlockInfo& info, BlockType& o) const {
    // FD coefficients for first and second derivative
    const Real h = info.h_gridpoint;
    const int thrID = omp_get_thread_num();
    const size_t nAgents = nAgentsPerBlock * nBlocks;
    size_t lastRealAgent = info.blockID;

    // Deal with the real agents first
    for (size_t k = 0; k < nAgentsPerBlock; k++) {
      const int ix = o.iAgentX[k], iy = o.iAgentY[k], iz = o.iAgentZ[k];
      const size_t agentID = nAgentsPerBlock*info.blockID + k;
      const std::vector<double> state = getState_uniform(lab, ix, iy, iz, h);
      if (&upcxx::current_persona()==&upcxx::backend::master)
        upcxx::liberate_master_persona();
        #pragma omp critical
        {
          upcxx::persona_scope scope(upcxx::master_persona());
          if (!timeOut){
            std::vector<double> Cs2_RL =
                comm.computeAction_upcxx(rlSeqInfo, state, reward, agentID);
            o(ix, iy, iz).chi = Cs2_RL[0];
          }
          else comm.truncSeq_upcxx(state, reward, agentID);
        }
    }

    // Then the fake agents
    for (int iz = 0; iz < FluidBlock::sizeZ; ++iz)
    for (int iy = 0; iy < FluidBlock::sizeY; ++iy)
    for (int ix = 0; ix < FluidBlock::sizeX; ++ix) {
      // one element per block is a proper agent: will add seq to train data
      // other are nThreads and are only there for thread safety
      // states get overwritten

      // LES coef can be stored in chi as long as we do not have obstacles
      // otherwise we will have to figure out smth
      //const bool bAgentTracked = evalStep and (ix == o.iAgentX && iy == o.iAgentY && iz == o.iAgentZ);
      const int localAgentID = getAgentId(ix, iy, iz, o.iAgentX, o.iAgentY, o.iAgentZ);

      //std::cout<<"Local Agent id"<< localAgentID << std::endl;
      if (localAgentID>=0) {
        lastRealAgent = nAgentsPerBlock*info.blockID + localAgentID;
        continue;
      }
      const size_t agentID =  nAgents + thrID;
      const std::vector<double> state = getState_uniform (lab, ix, iy, iz, h);

      if (!timeOut){
        std::vector<double> Cs2_RL =
            comm.computeAction_fakeAgent(rlSeqInfo, state, reward, agentID, lastRealAgent);
        o(ix, iy, iz).chi = Cs2_RL[0];
      }
      else comm.truncSeq_upcxx(state, reward, agentID);
    }
  }
};

template <typename Kernel>
void SGS_RL::runKernel(const Kernel& kernel)
{
  cubism::SynchronizerMPI<Real>& Synch = grid->sync(kernel);
  const int nthreads = omp_get_max_threads();
  LabMPI * labs = new LabMPI[nthreads];
  #pragma omp parallel for schedule(static, 1)
  for(int i = 0; i < nthreads; ++i) {
    labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
    labs[i].prepare(* sim.grid, Synch);
  }

  MPI_Barrier(grid->getCartComm());
  std::vector<cubism::BlockInfo> avail0 = Synch.avail_inner();
  std::vector<cubism::BlockInfo> avail1 = Synch.avail_halo();


  for(int i=0; i<sim.nprocs; ++i)
  {
    if(i == sim.rank)
    {
      int tid = 0;
      LabMPI& lab = labs[tid];
      const cubism::BlockInfo& I = avail0.size()? avail0[0] : avail1[0];
      FluidBlock& b = *(FluidBlock*)I.ptrBlock;
      lab.load(I, 0);
      kernel(lab, I, b);
    }
    MPI_Barrier(grid->getCartComm());
  }

  const int startInner = avail0.size()? 1 : 0;
  const int Ninner = avail0.size();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(static)
    for(int i=startInner; i<Ninner; i++) {
      const cubism::BlockInfo I = avail0[i];
      FluidBlock& b = *(FluidBlock*)I.ptrBlock;
      lab.load(I, 0);
      kernel(lab, I, b);
    }
  }

  const int Nhalo = avail1.size();
  const int startHalo = avail0.size()? 0 : 1;

  #pragma omp parallel
  {
  int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(static)
    for(int i=startHalo; i<Nhalo; i++) {
      const cubism::BlockInfo I = avail1[i];
      FluidBlock& b = *(FluidBlock*)I.ptrBlock;
      lab.load(I, 0);
      kernel(lab, I, b);
    }
  }

  if(labs != nullptr) {
    delete [] labs;
    labs = nullptr;
  }

  MPI_Barrier(grid->getCartComm());
}

SGS_RL::SGS_RL(SimulationData& s, Communicator* _comm, const int _step,
               const bool _timeOut, const bool _evalStep, const double _reward,
               const int _nAgentsPerBlock)
    : Operator(s), comm(_comm), step(_step), timeOut(_timeOut), evalStep(_evalStep), reward(_reward), nAgentsPerBlock(_nAgentsPerBlock) {}

void SGS_RL::operator()(const double dt) {
  sim.startProfiler("SGS_RL");
  const KernelSGS_RL K_SGS_RL(*comm, step, timeOut, evalStep, reward, sim.vInfo().size(), nAgentsPerBlock);

  runKernel<KernelSGS_RL>(K_SGS_RL);
  sim.stopProfiler();
  check("SGS_RL");
}

CubismUP_3D_NAMESPACE_END
