//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#ifdef _ACCFFT_
#include "SpectralManipACC.h"
#endif
#include "SpectralManipFFTW.h"

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

#include <sstream>
#include <iomanip>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

void HITstatistics::getTargetSpectrumFit(const Real _eps, const Real _nu,
                          std::vector<Real>& K, std::vector<Real>& E) const
{
  assert(K.size() == E.size() && K.size() > 0);
  // const Real gradFit = 0.8879967 std::sqrt(_eps) / std::sqrt(_nu);
  const Real LintFit = getIntegralLengthFit(_eps, _nu);
  const Real Lkolmogorov = getKolmogorovL(_eps, _nu);
  const Real C  = 5.7;
  const Real CI = 1e-3;
  const Real CE = 0;
  const Real BE = 4.9; // 5.2 from theory
  const Real P0 = 600; // should be 2, but we force large scales

  for (size_t i=0; i<K.size(); ++i) {
    K[i] = (i+0.5) * 2 * M_PI / L;
    const Real KI = K[i] * LintFit, KE4 = std::pow(K[i] * Lkolmogorov,4);
    const Real FL = std::pow(KI / (KI + CI), 5/3.0 + P0 );
    const Real FE = std::exp(-BE*(std::pow(KE4 + std::pow(CE,4), 0.25 ) -CE));
    E[i] = C * std::pow(_eps, 2/3.0) * std::pow(K[i], -5/3.0) * FL * FE;
  }
}

Real HITstatistics::getSimpleSpectrumFit(const Real _k, const Real _eps) const
{
  const Real C  = 5.7;
  return C * std::pow(_eps, 2/3.0) * std::pow(_k, -5/3.0);
}

Real HITstatistics::getIntegralTimeFit(const Real eps, const Real nu)
{
  return 0.93931475 * std::pow(eps, -1/3.0) * std::pow(nu, 1/6.0);
}

Real HITstatistics::getIntegralLengthFit(const Real eps, const Real nu)
{
  return 0.74885397 * std::pow(eps, -0.0233311) * std::pow(nu, 0.07192009);
}

Real HITstatistics::getTaylorMicroscaleFit(const Real eps, const Real nu)
{
  return 5.35507603 * std::pow(eps, -1/6.0) * std::sqrt(nu);
}

Real HITstatistics::getTurbKinEnFit(const Real eps, const Real nu)
{
  return 2.81574396 * std::pow(eps, 2/3.0);
}

Real HITstatistics::getHITReynoldsFit(const Real eps, const Real nu)
{
  return 7.33972668 * std::pow(eps, 1/6.0) / std::sqrt(nu);
}

Real HITstatistics::getKolmogorovL(const Real eps, const Real nu)
{
  return std::pow(eps, -0.25) * std::pow(nu, 0.75);
}

Real HITstatistics::getKolmogorovT(const Real eps, const Real nu)
{
  return std::sqrt(nu / eps);
}

EnergySpectrum::EnergySpectrum(const std::vector<Real>& _k,
                               const std::vector<Real>& _E) : k(_k), E(_E) {}
EnergySpectrum::EnergySpectrum(const std::vector<Real>& _k,
                               const std::vector<Real>& _E,
                               const std::vector<Real>& _s2) :
                               k(_k), E(_E), sigma2(_s2) {}

Real EnergySpectrum::interpE(const Real _k) const
{
  int idx_k = -1;
  int size = k.size();
  Real energy = 0.;
  for (int i = 0; i < size; ++i){
    if ( _k < k[i]){
      idx_k = i;
      break;
    }
  }
  if (idx_k > 0)
    energy = E[idx_k] + (E[idx_k-1] - E[idx_k]) / (k[idx_k] - k[idx_k-1]) * (k[idx_k] - _k);
  else if (idx_k==0)
    energy = E[idx_k] - (E[idx_k+1] - E[idx_k]) / (k[idx_k+1] - k[idx_k]) * (k[idx_k+1] - _k);
  else // idx_k==-1
    energy = E[size-1] + (E[size-2] - E[size-1]) / (k[size-1] - k[size-2]) * (k[size-1] - _k);

  energy = (energy < 0) ? 0. : energy;
  return energy;
}

Real EnergySpectrum::interpSigma2(const Real _k) const
{
  int idx_k = -1;
  int size = k.size();
  Real s2 = 0.;
  for (int i = 0; i < size; ++i){
    if ( _k < k[i]){
      idx_k = i;
      break;
    }
  }
  if (idx_k > 0)
    s2 = sigma2[idx_k] + (sigma2[idx_k-1] - sigma2[idx_k]) / (k[idx_k] - k[idx_k-1]) * (k[idx_k] - _k);
  else if (idx_k==0)
    s2 = sigma2[idx_k] - (sigma2[idx_k+1] - sigma2[idx_k]) / (k[idx_k+1] - k[idx_k]) * (k[idx_k] - _k);
  else // idx_k==-1
    s2 = sigma2[size-1] + (sigma2[size-2] - sigma2[size-1]) / (k[size-1] - k[size-2]) * (k[size-1] - _k);

  s2 = (s2 < 0) ? 0. : s2;
  return s2;
}

void EnergySpectrum::dump2File(const int nBin, const int nGrid, const Real lBox)
{
  const int nBins = std::ceil(std::sqrt(3)*nGrid)/2.0;//*waveFact);
  const Real binSize = M_PI*std::sqrt(3)*nGrid/(nBins*lBox);
  std::stringstream ssR;
  ssR<<"initialSpectrum.dat";
  std::ofstream f;
  f.open (ssR.str());
  f << std::left << std::setw(20) << "k" << std::setw(20) << "Pk" <<std::endl;
  for (int i = 0; i < nBins; ++i){
    const Real k_msr = (i+0.5)*binSize;
    f << std::left << std::setw(20) << k_msr << std::setw(20) << interpE(k_msr) <<std::endl;
  }
}

void initSpectralAnalysisSolver(SimulationData & sim)
{
  if(sim.spectralManip not_eq nullptr) return;
  if(not sim.bUseFourierBC) {
    printf("ERROR: spectral analysis functions support all-periodic BCs!\n");
    fflush(0); MPI_Abort(sim.app_comm, 1);
  }
  //#ifdef _ACCFFT_
    sim.spectralManip = new SpectralManipACC(sim);
  //#else
  //  sim.spectralManip = new SpectralManipFFTW(sim);
  //#endif
}

SpectralManip* initFFTWSpectralAnalysisSolver(SimulationData & sim)
{
  if(not sim.bUseFourierBC) {
    printf("ERROR: spectral analysis functions support all-periodic BCs!\n");
    fflush(0); MPI_Abort(sim.app_comm, 1);
  }
  return new SpectralManipFFTW(sim);
}

SpectralManip::SpectralManip(SimulationData & s) : sim(s)
{
  printf("New SpectralManip\n");
  int supported_threads;
  MPI_Query_thread(&supported_threads);
  if (supported_threads<MPI_THREAD_FUNNELED) {
    fprintf(stderr, "SpectralManip ERROR: MPI implementation does not support threads.\n");
    fflush(0); exit(1);
  }
}

SpectralManip::~SpectralManip()
{
}

CubismUP_3D_NAMESPACE_END
#undef MPIREAL
