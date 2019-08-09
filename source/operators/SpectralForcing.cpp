//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "SpectralForcing.h"
#include "SpectralManip.h"
#include "SpectralAnalysis.h"
#include "../utils/BufferedLogger.h"

CubismUP_3D_NAMESPACE_BEGIN using namespace cubism;

class KernelSpectralForcing
{
  typedef typename FluidGridMPI::BlockType BlockType;

 protected:
  const Real dt;
  const SpectralManip& sM;

 public:
  Real eps = 0, tke = 0, tkeFiltered = 0;

  void _cub2fftw() const;
  void _compute();
  void _fftw2cub(const Real factor) const;

  KernelSpectralForcing(SimulationData & s);
  ~KernelSpectralForcing() {}
};

KernelSpectralForcing::KernelSpectralForcing(SimulationData & s)
  : dt(s.dt), sM(*s.spectralManip)
{}

void KernelSpectralForcing::_cub2fftw() const
{
  const size_t NlocBlocks = sM.local_infos.size();
  Real * const data_u = sM.data_u;
  Real * const data_v = sM.data_v;
  Real * const data_w = sM.data_w;

  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM.local_infos[i].ptrBlock;
    const size_t offset = sM._offset( sM.local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix) {
      const size_t src_index = sM._dest(offset, iz, iy, ix);
      data_u[src_index] = b(ix,iy,iz).u;
      data_v[src_index] = b(ix,iy,iz).v;
      data_w[src_index] = b(ix,iy,iz).w;
    }
  }
}

void KernelSpectralForcing::_compute()
{
  fft_c *const cplxData_u = (fft_c *) sM.data_u;
  fft_c *const cplxData_v = (fft_c *) sM.data_v;
  fft_c *const cplxData_w = (fft_c *) sM.data_w;
  const Real waveFactorX = sM.waveFactor[0];
  const Real waveFactorY = sM.waveFactor[1];
  const Real waveFactorZ = sM.waveFactor[2];
  const long nKx = sM.nKx, nKy = sM.nKy, nKz = sM.nKz;
  const long sizeX = sM.gsize[0], sizeZ_hat = sM.nz_hat;
  const long loc_n1 = sM.local_n1, shifty = sM.shifty;

  #pragma omp parallel for reduction(+: eps, tke, tkeFiltered) schedule(static)
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
    const Real k2 = kx*kx + ky*ky + kz*kz, k_norm = std::sqrt(k2);

    const Real mult = (k==0) or (k==nKz/2) ? 1 : 2;
    const Real E = mult/2 * ( pow2_cplx(cplxData_u[linidx])
      + pow2_cplx(cplxData_v[linidx]) + pow2_cplx(cplxData_w[linidx]) );
    tke += E; // Total kinetic energy
    eps += k2 * E; // Dissipation rate
    if (k_norm > 0 && k_norm <= 2) {
      tkeFiltered += E;
    } else {
      cplxData_u[linidx][0] = 0;
      cplxData_u[linidx][1] = 0;
      cplxData_v[linidx][0] = 0;
      cplxData_v[linidx][1] = 0;
      cplxData_w[linidx][0] = 0;
      cplxData_w[linidx][1] = 0;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &tkeFiltered, 1, MPIREAL, MPI_SUM, sM.m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, sM.m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPIREAL, MPI_SUM, sM.m_comm);

  tkeFiltered =  tkeFiltered / pow2(sM.normalizeFFT);
  tke = tke / pow2(sM.normalizeFFT);
  eps = eps * 2*(sM.sim.nu) / pow2(sM.normalizeFFT);
}

void KernelSpectralForcing::_fftw2cub(const Real factor) const
{
  const size_t NlocBlocks = sM.local_infos.size();
  const Real * const data_u = sM.data_u;
  const Real * const data_v = sM.data_v;
  const Real * const data_w = sM.data_w;

  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM.local_infos[i].ptrBlock;
    const size_t offset = sM._offset( sM.local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; ++iz)
    for(int iy=0; iy<BlockType::sizeY; ++iy)
    for(int ix=0; ix<BlockType::sizeX; ++ix) {
      const size_t src_index = sM._dest(offset, iz, iy, ix);
      b(ix,iy,iz).u += factor * data_u[src_index];
      b(ix,iy,iz).v += factor * data_v[src_index];
      b(ix,iy,iz).w += factor * data_w[src_index];
    }
  }
}

SpectralForcing::SpectralForcing(SimulationData & s) : Operator(s)
{
  if (s.spectralManip == nullptr) s.spectralManip = new SpectralManip(s);
  sim.spectralManip->prepareFwd();
  sim.spectralManip->prepareBwd();
}

void SpectralForcing::operator()(const double dt)
{
  sim.startProfiler("SpectralForcing");
  SpectralManip * const sM = sim.spectralManip;
  assert(sM not_eq nullptr);
  KernelSpectralForcing K(sim);

  K._cub2fftw();

  sM->runFwd();

  K._compute();

  sM->runBwd();

  totalKinEn = K.tke;
  viscousDissip = K.eps;
  largeModesKinEn = K.tkeFiltered;
  sim.dissipationRate = (totalKinEnPrev - totalKinEn) / dt;
  totalKinEnPrev = totalKinEn;
  Real injectionRate = 0;
  //With non spectral IC, the target tke may not be defined here
  if      (sim.turbKinEn_target > 0) // inject energy to match target tke
       injectionRate = (sim.turbKinEn_target - totalKinEn)/dt;
  else if (sim.enInjectionRate  > 0) // constant power input:
       injectionRate =  sim.enInjectionRate;

  // If there's too much energy, let dissipation do its job
  if(sim.rank == 0) {
    std::stringstream &ssF = logger.get_stream("forcingData.dat");
    const std::string tab("\t");
    if(sim.step==0) {
      ssF<<"step \t time \t dt \t totalKinEn \t largeModesKinEn \t "\
           "viscousDissip \t totalDissipRate \t injectionRate\n";
    }

    ssF << sim.step << tab;
    ssF.setf(std::ios::scientific);
    ssF.precision(std::numeric_limits<float>::digits10 + 1);
    ssF<<sim.time<<tab<<sim.dt<<tab<<totalKinEn<<tab<<largeModesKinEn<<tab
       <<viscousDissip<<tab<<sim.dissipationRate<<tab<<injectionRate<<"\n";
  }

  const Real fac = sim.dt *injectionRate/(2*largeModesKinEn)/sM->normalizeFFT;

  if(fac>0) K._fftw2cub(fac);

  sim.stopProfiler();

  check("SpectralForcing");
}


CubismUP_3D_NAMESPACE_END
#undef MPIREAL
