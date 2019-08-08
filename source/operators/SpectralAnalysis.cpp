//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "operators/SpectralAnalysis.h"
#include "operators/SpectralManip.h"

#include <sys/stat.h>
#include <iomanip>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

SpectralAnalysis::SpectralAnalysis(SimulationData & s)
{
  if (s.spectralManip == nullptr) s.spectralManip = new SpectralManip(s);

  sM = s.spectralManip;

  nyquist = (int) sM->maxGridSize/2;
  nBin = nyquist-1;

  k_msr  = (Real*) std::calloc(nBin, sizeof(Real));
  for (int i = 0; i<nBin; i++) k_msr[i] = (i+1) * 2*M_PI / sM->maxBoxLength;

  E_msr  = (Real*) std::calloc(nBin, sizeof(Real));

  bComputeCs2Spectrum = s.bComputeCs2Spectrum;
  if (bComputeCs2Spectrum)
    cs2_msr  = (Real*) std::calloc(nBin, sizeof(Real));

  sM->prepareFwd();
}

void SpectralAnalysis::_cub2fftw()
{
  // Let's also compute u_avg here
  const size_t NlocBlocks = sM->local_infos.size();
  Real * const data_u = sM->data_u;
  Real * const data_v = sM->data_v;
  Real * const data_w = sM->data_w;
  Real * const data_cs2 = sM->data_cs2;
  const SpectralManip & helper = * sM;

  u_avg[0] = 0; u_avg[1] = 0; u_avg[2] = 0; unorm_avg = 0;
  #pragma omp parallel for reduction(+: u_avg[:3], unorm_avg) schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i)
  {
    const BlockType& b = *(BlockType*) helper.local_infos[i].ptrBlock;
    const size_t offset = helper._offset( helper.local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix)
    {
      const size_t src_index = helper._dest(offset, iz, iy, ix);
      const Real u = b(ix,iy,iz).u;
      const Real v = b(ix,iy,iz).v;
      const Real w = b(ix,iy,iz).w;
      u_avg[0] += u; u_avg[1] += v; u_avg[2] += w;
      unorm_avg += (u*u + v*v + w*w)/2;
      data_u[src_index] = u;
      data_v[src_index] = v;
      data_w[src_index] = w;
      if (bComputeCs2Spectrum) data_cs2[src_index] = b(ix,iy,iz).chi;
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
  fft_c *const cplxData_cs2  = (fft_c *) sM->data_cs2;
  const Real waveFactorX = sM->waveFactor[0];
  const Real waveFactorY = sM->waveFactor[1];
  const Real waveFactorZ = sM->waveFactor[2];
  const long nKx = sM->nKx, nKy = sM->nKy, nKz = sM->nKz;
  const long loc_n1 = sM->local_n1, shifty = sM->shifty;
  const long sizeX = sM->gsize[0], sizeZ_hat = sM->nz_hat;

  tke = 0;
  eps = 0;
  tau_integral = 0;
  memset(E_msr, 0, nBin * sizeof(Real));

  // Let's only measure spectrum up to Nyquist.
  #pragma omp parallel for reduction(+ : E_msr[:nBin], tke, eps, tau_integral)  schedule(static)
  for(long j = 0; j<loc_n1; ++j)
  for(long i = 0; i<sizeX;  ++i)
  for(long k = 0; k<sizeZ_hat; ++k)
  {
    const long linidx = (j*sizeX +i)*sizeZ_hat + k;
    const long ii = (i <= nKx/2) ? i : -(nKx-i);
    const long l = shifty + j; //memory index plus shift due to decomp
    const long jj = (l <= nKy/2) ? l : -(nKy-l);
    const long kk = (k <= nKz/2) ? k : -(nKz-k);

    const Real kx = ii*waveFactorX, ky = jj*waveFactorY, kz = kk*waveFactorZ;
    const Real mult = (k==0) or (k==nKz/2) ? 1 : 2;
    const Real k2 = kx*kx + ky*ky + kz*kz, k_norm = std::sqrt(k2);
    const Real ks = std::sqrt(ii*ii + jj*jj + kk*kk);
    const Real E = mult/2 * ( pow2_cplx(cplxData_u[linidx])
      + pow2_cplx(cplxData_v[linidx]) + pow2_cplx(cplxData_w[linidx]) );

    // Total kinetic energy
    tke += E;
    // Dissipation rate
    eps += k2 * E;
    // Large eddy turnover time
    tau_integral += (k_norm > 0) ? E / k_norm : 0.;
    if (ks <= nyquist) {
      int binID = std::floor(ks * (nyquist-1)/nyquist);
      E_msr[binID] += E;
      if (bComputeCs2Spectrum){
        const Real cs2 = std::sqrt(pow2_cplx(cplxData_cs2[linidx]));
        cs2_msr[binID] += mult*cs2;
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, E_msr,  nBin, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tau_integral, 1, MPIREAL, MPI_SUM, sM->m_comm);

  if (bComputeCs2Spectrum){
    assert(false);
    //#pragma omp parallel reduction(+ : cs2_msr[:nBin])
    //MPI_Allreduce(MPI_IN_PLACE, cs2_msr, nBin, MPIREAL, MPI_SUM, sM->m_comm);
  }

  const size_t normalize = pow2(sM->normalizeFFT);

  for (int binID = 0; binID < nBin; binID++) {
    E_msr[binID] /= normalize;
    if (bComputeCs2Spectrum)
      cs2_msr[binID] /= sM->normalizeFFT;
  }

  tke = tke / normalize;
  eps = eps * 2*(sM->sim.nu)/ normalize;

  uprime = std::sqrt(2*tke/3);
  lambda = std::sqrt(15.0*sM->sim.nu/eps)*uprime;
  Re_lambda = uprime*lambda/sM->sim.nu;

  tau_integral = M_PI / (2*pow3(uprime)) * tau_integral / normalize;
}

void SpectralAnalysis::run()
{
  _cub2fftw();
  sM->runFwd();
  _compute();
}

void SpectralAnalysis::dump2File(const int nFile) const
{
  std::stringstream ssR;
  ssR << "analysis/spectralAnalysis_" << std::setfill('0') << std::setw(9)
      << nFile;
  std::ofstream f;
  f.open(ssR.str());
  f << std::left << "Spectral Analysis :" << "\n";
  f << std::left << std::setw(15) << "time"
    << std::setw(15) << sM->sim.time
    << " #simulation time" << "\n";

  f << std::left << std::setw(15) << "lBox"
    << std::setw(15) << sM->maxBoxLength
    << " #simulation box length" << "\n";

  f << std::left << std::setw(15) << "tke"
    << std::setw(15) << tke
    << " #turbulent kinetic energy" << "\n";

  f << std::left << std::setw(15) << "eps"
    << std::setw(15) << eps
    << " #Viscous dissipation rate" << "\n";

  f << std::left << std::setw(15) << "eps_f"
    << std::setw(15) << sM->sim.dissipationRate
    << " #Total dissipation rate" << "\n";

  f << std::left << std::setw(15) << "lambda"
    << std::setw(15) << lambda
    << " #Taylor microscale" << "\n";

  f << std::left << std::setw(15) << "Re_lambda"
    << std::setw(15) << Re_lambda
    << " #Turbulent Reynolds based on lambda" << "\n";

  f << std::left << std::setw(15) << "eta"
    << std::setw(15) << std::pow(pow3(sM->sim.nu)/eps, 1.0/4)
    << " #Kolmogorov length scale" << "\n";

  f << std::left << std::setw(15) << "tau_integral"
    << std::setw(15) << tau_integral
    << " #Integral time scale" << "\n";

  f << std::left << std::setw(15) << "tau_eta"
    << std::setw(15) << std::sqrt(sM->sim.nu/eps)
    << " #Kolmogorov time scale" << "\n";

  f << std::left << std::setw(15) << "nu_sgs"
    << std::setw(15) << sM->sim.nu_sgs
    << " #Average SGS viscosity" << "\n";

  f << std::left << std::setw(15) << "cs2_avg"
    << std::setw(15) << sM->sim.cs2_avg
    << " #Average Cs2 if dynamic model" << "\n";

  f << std::left << std::setw(15) << "mean_grad"
    << std::setw(15) << sM->sim.grad_mean
    << " #Average gradient magnitude" << "\n";

  f << std::left << std::setw(15) << "std_grad"
    << std::setw(15) << sM->sim.grad_std
    << " #Stdev gradient magnitude" << "\n"
    << "\n";

  f << std::left << std::setw(15) << "k * (lBox/2pi)"
    << std::setw(15) << "E_k" << "\n";

  for (int i = 0; i < nBin; i++)
    f << std::left << std::setw(15) << k_msr[i]
                   << std::setw(15) << E_msr[i] << "\n";
  f.flush();
  f.close();

  if (bComputeCs2Spectrum)
  {
    std::stringstream ssR_cs2;
    ssR_cs2 << "analysis/spectralAnalysisCs2_" << std::setfill('0')
            << std::setw(9) << nFile;
    f.open(ssR_cs2.str());
    f << std::left << "Cs2 spectrum :" << "\n";
    f << std::left << std::setw(15) << "k * (lBox/2pi)"
                   << std::setw(15) << "Cs2_k" << "\n";
    for (int i = 0; i < nBin; i++)
      f << std::left << std::setw(15) << k_msr[i]
                     << std::setw(15) << cs2_msr[i] << "\n";
    f.flush();
    f.close();
  }
}

void SpectralAnalysis::reset()
{
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
