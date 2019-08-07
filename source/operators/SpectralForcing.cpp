//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "operators/SpectralForcing.h"
#include "operators/SpectralManip.h"
#include "operators/SpectralAnalysis.h"

CubismUP_3D_NAMESPACE_BEGIN using namespace cubism;

class KernelSpectralForcing
{
  typedef typename FluidGridMPI::BlockType BlockType;
 protected:
  const Real dt;
 public:
  Real tke_filtered = 0.0;
  Real tke = 0.0;

  void _cub2fftw() const;
  void _compute();
  void _fftw2cub();

  SpectralManip *sM;

  KernelSpectralForcing(SimulationData & s);
  ~KernelSpectralForcing() {}
  void run();
};

KernelSpectralForcing::KernelSpectralForcing(SimulationData & s) : dt(s.dt)
{
  if (s.spectralManip == nullptr)
    s.spectralManip = new SpectralManip(s);

  sM = s.spectralManip;
  sM->prepareFwd();
  sM->prepareBwd();
}

void KernelSpectralForcing::_cub2fftw() const
{
  const size_t NlocBlocks = sM->local_infos.size();
  Real * const data_u = sM->data_u;
  Real * const data_v = sM->data_v;
  Real * const data_w = sM->data_w;

  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      data_u[src_index] = b(ix,iy,iz).u;
      data_v[src_index] = b(ix,iy,iz).v;
      data_w[src_index] = b(ix,iy,iz).w;
    }
  }
}

void KernelSpectralForcing::_compute()
{
  fft_c *const cplxData_u = (fft_c *) sM->data_u;
  fft_c *const cplxData_v = (fft_c *) sM->data_v;
  fft_c *const cplxData_w = (fft_c *) sM->data_w;
  const Real waveFactorX = sM->waveFactor[0];
  const Real waveFactorY = sM->waveFactor[1];
  const Real waveFactorZ = sM->waveFactor[2];
  const long nKx = sM->nKx, nKy = sM->nKy, nKz = sM->nKz;
  const long sizeX = sM->gsize[0], sizeZ_hat = sM->nz_hat;
  const long loc_n1 = sM->local_n1, shifty = sM->shifty;

  #pragma omp parallel for reduction(+: tke, tke_filtered) schedule(static)
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
    const Real k2 = kx*kx + ky*ky + kz*kz;
    const Real k_norm = std::sqrt(k2);

    const int mult = (k==0) or (k==nKz/2) ? 1 : 2;
    const Real E = pow2_cplx(cplxData_u[linidx])
                 + pow2_cplx(cplxData_v[linidx])
                 + pow2_cplx(cplxData_w[linidx]);
    tke += mult*E/2;
    if (k_norm > 0 && k_norm <= 2) {
      tke_filtered += mult*E/2;
    } else {
      cplxData_u[linidx][0] = 0;
      cplxData_u[linidx][1] = 0;
      cplxData_v[linidx][0] = 0;
      cplxData_v[linidx][1] = 0;
      cplxData_w[linidx][0] = 0;
      cplxData_w[linidx][1] = 0;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &tke_filtered, 1, MPIREAL, MPI_SUM, sM->m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, sM->m_comm);

  tke_filtered /=  pow2(sM->normalizeFFT);
  tke /=  pow2(sM->normalizeFFT);
}

void KernelSpectralForcing::_fftw2cub()
{
  //With non spectral IC, the target tke may not be defined here
  if      (sM->sim.turbKinEn_target > 0) // inject energy to match target tke
       sM->sim.injectedPower = (sM->sim.turbKinEn_target - tke)/dt;
  else if (sM->sim.enInjectionRate > 0) // constant power input:
       sM->sim.injectedPower =  sM->sim.enInjectionRate;
  else sM->sim.injectedPower = 0;

  // If there's too much energy, let dissipation do its job
  if (sM->sim.injectedPower <= 0) {
    sM->sim.injectedPower = 0;
    return;
  }

  printf("Total Kin E = %f, E_|k|<2 = %f, injected power %f (energy = %f)\n",
    tke, tke_filtered, sM->sim.injectedPower, dt * sM->sim.injectedPower);
  const Real F = dt * sM->sim.injectedPower /(2*tke_filtered) /sM->normalizeFFT;
  const size_t NlocBlocks = sM->local_infos.size();
  const Real * const data_u = sM->data_u;
  const Real * const data_v = sM->data_v;
  const Real * const data_w = sM->data_w;

  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      b(ix,iy,iz).u += F * data_u[src_index];
      b(ix,iy,iz).v += F * data_v[src_index];
      b(ix,iy,iz).w += F * data_w[src_index];
    }
  }
}

void KernelSpectralForcing::run()
{
  _cub2fftw();
  sM->runFwd();
  _compute();
  sM->runBwd();
  _fftw2cub();
}

void SpectralForcing::operator()(const double dt)
{
  sim.startProfiler("SpectralForcing");
  KernelSpectralForcing K(sim);
  K.run();
  sim.stopProfiler();
  check("SpectralForcing");
}


CubismUP_3D_NAMESPACE_END
#undef MPIREAL
