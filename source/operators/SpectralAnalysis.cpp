//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "operators/SpectralAnalysis.h"

#include <sys/stat.h>
#include <iomanip>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

SpectralAnalysis::SpectralAnalysis(SimulationData & s)
{
  if (s.spectralManip == nullptr)
    s.spectralManip = new SpectralManip(s);

  sM = s.spectralManip;

  nyquist = (int) sM->maxGridSize/2;
  nBin = nyquist-1;

  k_msr  = (Real*) std::calloc(nBin, sizeof(Real));
  for (int i = 0; i<nBin; i++)k_msr[i] = (i+1) * 2*M_PI / sM->maxBoxLength;

  E_msr  = (Real*) std::calloc(nBin, sizeof(Real));

  sM->prepareFwd();
}

void SpectralAnalysis::_cub2fftw()
{
  // Let's also compute u_avg here
  const size_t NlocBlocks = sM->local_infos.size();
  #pragma omp parallel for reduction(+: u_avg[:3])
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      const Real u = b(ix,iy,iz).u;
      const Real v = b(ix,iy,iz).v;
      const Real w = b(ix,iy,iz).w;
      u_avg[0] += u; u_avg[1] += v; u_avg[2] += w;
      sM->data_u[src_index] = u;
      sM->data_v[src_index] = v;
      sM->data_w[src_index] = w;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, u_avg, 3, MPIREAL, MPI_SUM, sM->m_comm);
  // normalizeFFT is the total number of grid cells
  u_avg[0] /= sM->normalizeFFT;
  u_avg[1] /= sM->normalizeFFT;
  u_avg[2] /= sM->normalizeFFT;
}

void SpectralAnalysis::_compute()
{
  fft_c *const cplxData_u  = (fft_c *) sM->data_u;
  fft_c *const cplxData_v  = (fft_c *) sM->data_v;
  fft_c *const cplxData_w  = (fft_c *) sM->data_w;
  // Let's only measure spectrum up to Nyquist.
#pragma omp parallel for reduction(+ : E_msr[:nBin], tke, eps, tau_integral)
  for(long j = 0; j<static_cast<long>(sM->local_n1); ++j)
  for(long i = 0; i<static_cast<long>(sM->gsize[0]); ++i)
  for(long k = 0; k<static_cast<long>(sM->nz_hat);   ++k) {
    const size_t linidx = (j*sM->gsize[0] +i)*sM->nz_hat + k;
    const long ii = (i <= sM->nKx/2) ? i : -(sM->nKx-i);
    const long l = sM->shifty + j; //memory index plus shift due to decomp
    const long jj = (l <= sM->nKy/2) ? l : -(sM->nKy-l);
    const long kk = (k <= sM->nKz/2) ? k : -(sM->nKz-k);

    const Real E = pow2_cplx(cplxData_u[linidx])
                 + pow2_cplx(cplxData_v[linidx])
                 + pow2_cplx(cplxData_w[linidx]);
    const Real kx = ii*sM->waveFactor[0], ky = jj*sM->waveFactor[1], kz = kk*sM->waveFactor[2];
    const Real k_norm = sqrt(kx*kx + ky*ky + kz*kz);
    const Real ks = sqrt(ii*ii + jj*jj + kk*kk);
    tke += E;
    eps += pow2(k_norm) * E;
    tau_integral += (k_norm > 0) ? E / k_norm : 0.;
    if (ks <= nyquist){
      int binID = std::floor(ks * (nyquist-1)/nyquist);
      E_msr[binID] += E;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, E_msr,  nBin, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tau_integral, 1, MPIREAL, MPI_SUM, sM->m_comm);
  const size_t normalize = pow2(sM->normalizeFFT);
  for (int binID = 0; binID < nBin; binID++) E_msr[binID] /= normalize;

  tke = tke / normalize;
  eps = eps * 2*(sM->sim.nu + sM->sim.nu_sgs) / normalize;
  uprime = sqrt(2*tke/3);
  lambda = sqrt(15.0*sM->sim.nu/eps)*uprime;
  Re_lambda = uprime*lambda/sM->sim.nu;

  tau_integral = M_PI / (2*pow3(uprime)) * tau_integral / normalize;
}

void SpectralAnalysis::run()
{
  _cub2fftw();
  sM->runFwd();
  _compute();
}

void SpectralAnalysis::dump2File(const int nFile) const {
  std::stringstream ssR;
  ssR << "analysis/spectralAnalysis_" << std::setfill('0') << std::setw(9)
      << nFile;
  std::ofstream f;
  f.open(ssR.str());
  f << std::left << "Spectral Analysis :" << std::endl;
  f << std::left << std::setw(15) << "time" << std::setw(15) << sM->sim.time
    << " #simulation time" << std::endl;

  f << std::left << std::setw(15) << "lBox" << std::setw(15) << sM->maxBoxLength
    << " #simulation box length" << std::endl;

  f << std::left << std::setw(15) << "tke" << std::setw(15) << tke
    << " #turbulent kinetic energy" << std::endl;

  f << std::left << std::setw(15) << "eps" << std::setw(15) << eps
    << " #dissipation rate" << std::endl;

  f << std::left << std::setw(15) << "uprime" << std::setw(15) << uprime
    << " #Average RMS velocity" << std::endl;

  f << std::left << std::setw(15) << "lambda" << std::setw(15) << lambda
    << " #Taylor microscale" << std::endl;

  f << std::left << std::setw(15) << "Re_lambda" << std::setw(15) << Re_lambda
    << " #Turbulent Reynolds based on lambda" << std::endl;

  f << std::left << std::setw(15) << "eta" << std::setw(15) << pow(pow3(sM->sim.nu)/eps, 1.0/4)
    << " #Kolmogorov length scale" << std::endl;

  f << std::left << std::setw(15) << "tau_integral" << std::setw(15) << tau_integral
    << " #Integral time scale" << std::endl;

  f << std::left << std::setw(15) << "tau_eta" << std::setw(15) << sqrt(sM->sim.nu/eps)
    << " #Kolmogorov time scale" << std::endl;

  f << std::left << std::setw(15) << "nu_sgs" << std::setw(15) << sM->sim.nu_sgs
    << " #Average SGS viscosity" << std::endl;

  f << std::left << std::setw(15) << "cs2_rl" << std::setw(15) << sM->sim.cs2_rl
    << " #Average Cs2 from RL" << std::endl
    << std::endl;

  f << std::left << std::setw(15) << "k * (lBox/2pi)" << std::setw(15) << "E_k" << std::endl;
  for (int i = 0; i < nBin; i++) {
    f << std::left << std::setw(15) << k_msr[i] << std::setw(15) << E_msr[i] << std::endl;
  }
  f.close();
}

void SpectralAnalysis::reset(){
  memset(E_msr, 0, nBin * sizeof(Real));
  memset(k_msr, 0, nBin * sizeof(Real));
}

SpectralAnalysis::~SpectralAnalysis()
{
  delete[] k_msr;
  delete[] E_msr;
}

CubismUP_3D_NAMESPACE_END
#undef MPIREAL
