//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifndef CubismUP_3D_SpectralManip_h
#define CubismUP_3D_SpectralManip_h

#include "../SimulationData.h"
#include "../poisson/PoissonSolver_common.h"

#include <Cubism/BlockInfo.h>

#include <cassert>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <vector>

CubismUP_3D_NAMESPACE_BEGIN

template <class T>
inline T pow2(const T val) {
    return val*val;
}

template <class T>
inline T pow3(const T val) {
    return val*val*val;
}

inline Real pow2_cplx(const fft_c cplx_val) {
  return pow2(cplx_val[0]) + pow2(cplx_val[1]);
}

struct energySpectrum
{
public:
  const std::vector<Real> k;
  const std::vector<Real> E;
  const std::vector<Real> sigma2;

  energySpectrum(std::vector<Real> _k, std::vector<Real> _E) : k(_k), E(_E) {}
  energySpectrum(std::vector<Real> _k, std::vector<Real> _E,
                 std::vector<Real> _sigma2)
      : k(_k), E(_E), sigma2(_sigma2) {}

  Real interpE(const Real k);
  Real interpSigma2(const Real k);
  void dump2File(const int nBin, const int nGrid, const Real h);
};

class SpectralManip
{
  friend class SpectralIcGenerator;
  friend class SpectralAnalysis;
  friend class KernelSpectralForcing;
 private:
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;
  // MPI related
  const MPI_Comm m_comm = grid.getCartComm();
  const int m_rank = sim.rank, m_size = sim.nprocs;

  static constexpr int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const std::vector<cubism::BlockInfo> local_infos = grid.getResidentBlocksInfo();

  const size_t mybpd[3] = {
      static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
  };
  const size_t gsize[3] = {
      static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]),
      static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]),
      static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])
  };
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

  size_t stridez = 0;
  size_t stridey = 0;
  size_t stridex = 0;
  size_t data_size = 0;

  long shifty=0;

  const long nKx = static_cast<long>(gsize[0]);
  const long nKy = static_cast<long>(gsize[1]);
  const long nKz = static_cast<long>(gsize[2]);

  const long maxGridSize = std::max({gsize[0], gsize[1], gsize[2]});
  const long minGridSize = std::max({gsize[0], gsize[1], gsize[2]});
  const std::array<Real, 3> lBox = sim.extent;
  const Real maxBoxLength = std::max({lBox[0], lBox[1], lBox[2]});
  const Real minBoxLength = std::min({lBox[0], lBox[1], lBox[2]});

  const size_t nz_hat = gsize[2]/2+1;
  const double h = sim.uniformH();
  const double waveFactor[3] = {2.0 * M_PI / lBox[0],
                                2.0 * M_PI / lBox[1],
                                2.0 * M_PI / lBox[2]};
  bool bAllocFwd = false, bAllocBwd = false;
  bool bAllocCs2 = false;

  const int nBin = std::ceil(std::sqrt(3.0)*maxGridSize/2.0)+1;
  const Real binSize = M_PI*std::sqrt(3.0)*maxGridSize/(nBin*maxBoxLength);

  Real * data_u, * data_v, * data_w, * data_cs2;
  void * fwd_u, * fwd_v, * fwd_w, * fwd_cs2;
  void * bwd_u, * bwd_v, * bwd_w;

public:

  const size_t normalizeFFT = gsize[0] * gsize[1] * gsize[2];

  SpectralManip(SimulationData & s);
  ~SpectralManip();

  void prepareFwd();
  void prepareBwd();

  inline size_t _offset(const int blockID) const
  {
    const cubism::BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }
  inline size_t _offset_ext(const cubism::BlockInfo &info) const
  {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
  }
  inline size_t _offset(const cubism::BlockInfo &info) const
  {
    assert(stridez>0);
    assert(stridey>0);
    assert(stridex>0);
    assert(data_size>0);
    const int myIstart[3] = {
      info.index[0]*bs[0],
      info.index[1]*bs[1],
      info.index[2]*bs[2]
    };
    return stridez*myIstart[2] + stridey*myIstart[1] + stridex*myIstart[0];
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    assert(stridez>0);
    assert(stridey>0);
    assert(stridex>0);
    assert(data_size>0);
    return offset + stridez*z + stridey*y + stridex*x;
  }

  void runFwd() const;
  void runBwd() const;

  void reset() const
  {
    memset(data_u, 0, data_size * sizeof(Real));
    memset(data_v, 0, data_size * sizeof(Real));
    memset(data_w, 0, data_size * sizeof(Real));
  }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralManip_h
