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
  #pragma omp parallel for
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      sM->data_u[src_index] = b(ix,iy,iz).u;
      sM->data_v[src_index] = b(ix,iy,iz).v;
      sM->data_w[src_index] = b(ix,iy,iz).w;
    }
  }
}

void KernelSpectralForcing::_compute()
{
  fft_c *const cplxData_u = (fft_c *) sM->data_u;
  fft_c *const cplxData_v = (fft_c *) sM->data_v;
  fft_c *const cplxData_w = (fft_c *) sM->data_w;

#pragma omp parallel for reduction(+: tke, tke_filtered)
  for(long j = 0; j<static_cast<long>(sM->local_n1); ++j)
  for(long i = 0; i<static_cast<long>(sM->gsize[0]); ++i)
  for(long k = 0; k<static_cast<long>(sM->nz_hat);   ++k) {
    const size_t linidx = (j*sM->gsize[0] +i)*sM->nz_hat + k;
    const long ii = (i <= sM->nKx/2) ? i : -(sM->nKx-i);
    const long l = sM->shifty + j; //memory index plus shift due to decomp
    const long jj = (l <= sM->nKy/2) ? l : -(sM->nKy-l);
    const long kk = (k <= sM->nKz/2) ? k : -(sM->nKz-k);

    const Real kx = ii*sM->waveFactor[0], ky = jj*sM->waveFactor[1], kz = kk*sM->waveFactor[2];
    const Real k2 = kx*kx + ky*ky + kz*kz;
    const Real k_norm = sqrt(k2);

    const int mult = (k==0) or (k==sM->nKz/2) ? 1 : 2;
    const Real E = pow2_cplx(cplxData_u[linidx])
                 + pow2_cplx(cplxData_v[linidx])
                 + pow2_cplx(cplxData_w[linidx]);
    tke += 0.5*mult*E;
    if (0 < k_norm and k_norm <= 2) tke_filtered += 0.5*mult*E;
    else
    {
      cplxData_u[linidx][0] = cplxData_u[linidx][1] = 0.;
      cplxData_v[linidx][0] = cplxData_v[linidx][1] = 0.;
      cplxData_w[linidx][0] = cplxData_w[linidx][1] = 0.;
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
  if (sM->sim.tkeTgt ==0) sM->sim.tkeTgt = tke;

  const Real eps = (sM->sim.tkeTgt - tke)/dt;

  // If there's too much energy, let dissipation do its job
  if (eps < 0) {
    sM->sim.epsForcing = 0;
    return;
  }

  // Otherwise, inject energy to match target tke
  sM->sim.epsForcing = eps;
  const Real fac = dt * eps / (2*tke_filtered) / sM->normalizeFFT;
  const size_t NlocBlocks = sM->local_infos.size();
  #pragma omp parallel for
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      b(ix,iy,iz).u += fac*sM->data_u[src_index];
      b(ix,iy,iz).v += fac*sM->data_v[src_index];
      b(ix,iy,iz).w += fac*sM->data_w[src_index];
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
