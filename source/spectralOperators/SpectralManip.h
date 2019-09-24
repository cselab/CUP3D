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
#include "Cubism/BlockInfo.h"

#include <vector>
#include <cassert>
#include <cstring>

CubismUP_3D_NAMESPACE_BEGIN

template <class T>
inline T pow2(const T val) {
    return val*val;
}

template <class T>
inline T pow3(const T val) {
    return val*val*val;
}

class SpectralManip;

void initSpectralAnalysisSolver(SimulationData & sim);
SpectralManip* initFFTWSpectralAnalysisSolver(SimulationData & sim);

struct EnergySpectrum
{
  const std::vector<Real> k;
  const std::vector<Real> E;
  const std::vector<Real> sigma2;

  EnergySpectrum(const std::vector<Real> &_k, const std::vector<Real> &_E);
  EnergySpectrum(const std::vector<Real> &_k, const std::vector<Real> &_E,
                 const std::vector<Real> &_sigma2);

  Real interpE(const Real k) const;
  Real interpSigma2(const Real k) const;
  void dump2File(const int nBin, const int nGrid, const Real h);
};

struct HITstatistics
{
  HITstatistics() = delete;

  HITstatistics(const int maxGridSize, const Real maxBoxLength):
    N(maxGridSize), L(maxBoxLength),
    k_msr(new Real[nBin]), E_msr(new Real[nBin]), cs2_msr(new Real[nBin])
  {
    //printf("maxGridSize %d %d %d\n", maxGridSize, N, nyquist);
    reset();
    for (int i = 0; i<nBin; ++i) k_msr[i] = (i+1) * 2*M_PI / L;
  }

  HITstatistics(const HITstatistics&c) : N(c.N), L(c.L),
    k_msr(new Real[nBin]), E_msr(new Real[nBin]), cs2_msr(new Real[nBin])
  {
    //printf("maxGridSize %d %d %d\n", c.N, N, nyquist);
    reset();
    for (int i = 0; i<nBin; ++i) k_msr[i] = (i+1) * 2*M_PI / L;
  }

  ~HITstatistics()
  {
    delete [] k_msr;
    delete [] E_msr;
    delete [] cs2_msr;
  }

  void reset()
  {
    tke = 0; eps = 0; tau_integral = 0; tke_filtered = 0;
    lambda = 0; uprime = 0; Re_lambda = 0;
    memset(E_msr, 0, nBin * sizeof(Real));
    memset(cs2_msr, 0, nBin * sizeof(Real));
  }

  Real getSimpleSpectrumFit(const Real _k, const Real _eps) const;
  void getTargetSpectrumFit(const Real eps, const Real nu,
                            std::vector<Real>& K, std::vector<Real>& E) const;
  static Real getIntegralTimeFit(const Real eps, const Real nu);
  static Real getIntegralLengthFit(const Real eps, const Real nu);
  static Real getTaylorMicroscaleFit(const Real eps, const Real nu);
  static Real getHITReynoldsFit(const Real eps, const Real nu);
  static Real getTurbKinEnFit(const Real eps, const Real nu);
  static Real getKolmogorovL(const Real eps, const Real nu);
  static Real getKolmogorovT(const Real eps, const Real nu);

  // Parameters of the histogram
  const int N, nyquist = N/2, nBin = nyquist-1;
  const Real L;

  // Output of the analysis
  Real tke = 0, eps = 0, tau_integral = 0, l_integral = 0, tke_filtered = 0;
  Real lambda = 0, uprime = 0, Re_lambda = 0;
  Real * const k_msr;
  Real * const E_msr;
  Real * const cs2_msr;
};

class SpectralManip
{
  friend class SpectralIcGenerator;
  friend class SpectralAnalysis;
  friend class SpectralForcing;

protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;

  size_t stridez = 0;
  size_t stridey = 0;
  size_t stridex = 0;
  size_t data_size = 0;

  // spectral manips are only supported for fully periodic flows
  bool bAllocFwd = false;
  bool bAllocBwd = false;

  Real * data_u, * data_v, * data_w;
  // * data_cs2;

public:
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
  const Real normalizeFFT = gsize[0] * gsize[1] * gsize[2];
  const size_t nz_hat = gsize[2]/2+1;
  const long maxGridN = std::max({gsize[0], gsize[1], gsize[2]});
  const Real maxGridL = std::max({sim.extent[0], sim.extent[1], sim.extent[2]});
  HITstatistics stats = HITstatistics(maxGridN, maxGridL);
  //const double h = sim.uniformH();

  SpectralManip(SimulationData & s);
  virtual ~SpectralManip();

  virtual void prepareFwd() = 0;
  virtual void prepareBwd() = 0;

  size_t _offset(const int blockID) const
  {
    const cubism::BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }

  size_t _offset_ext(const cubism::BlockInfo &info) const
  {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
  }

  size_t _offset(const cubism::BlockInfo &info) const
  {
    assert(stridez>0);
    assert(stridey>0);
    assert(stridex>0);
    assert(data_size>0);
    const int myIstart[3] = {
      info.index[0] * bs[0],
      info.index[1] * bs[1],
      info.index[2] * bs[2]
    };
    return stridez*myIstart[2] + stridey*myIstart[1] + stridex*myIstart[0];
  }

  size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    assert(stridez>0);
    assert(stridey>0);
    assert(stridex>0);
    assert(data_size>0);
    return offset + stridez*z + stridey*y + stridex*x;
  }

  virtual void runFwd() const = 0;
  virtual void runBwd() const = 0;

  virtual void _compute_largeModesForcing() = 0;
  virtual void _compute_analysis() = 0;
  virtual void _compute_IC(const std::vector<Real> &K,
                           const std::vector<Real> &E) = 0;

  void reset() const
  {
    memset(data_u, 0, data_size * sizeof(Real));
    memset(data_v, 0, data_size * sizeof(Real));
    memset(data_w, 0, data_size * sizeof(Real));
  }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_SpectralManip_h
