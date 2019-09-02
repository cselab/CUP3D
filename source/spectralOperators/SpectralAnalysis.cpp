//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "SpectralAnalysis.h"
#include "SpectralManip.h"
#include <Cubism/HDF5Dumper_MPI.h>

#include <sys/stat.h>
#include <iomanip>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

SpectralAnalysis::SpectralAnalysis(SimulationData & s)
{
  initSpectralAnalysisSolver(s);
  s.spectralManip->prepareFwd();
  sM = s.spectralManip;
}

void SpectralAnalysis::_cub2fftw()
{
  // Let's also compute u_avg here
  const size_t NlocBlocks = sM->local_infos.size();
  Real * const data_u = sM->data_u;
  Real * const data_v = sM->data_v;
  Real * const data_w = sM->data_w;
  //Real * const data_cs2 = sM->data_cs2;
  const SpectralManip & helper = * sM;

  u_avg[0] = 0; u_avg[1] = 0; u_avg[2] = 0; unorm_avg = 0;
  #pragma omp parallel for reduction(+: u_avg[:3]) schedule(static)
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
      //unorm_avg += (u*u + v*v + w*w)/2;
      data_u[src_index] = u;
      data_v[src_index] = v;
      data_w[src_index] = w;
      // data_cs2[src_index] = b(ix,iy,iz).chi;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, u_avg, 3, MPIREAL, MPI_SUM, sM->m_comm);
  // normalizeFFT is the total number of grid cells
  u_avg[0] /= sM->normalizeFFT;
  u_avg[1] /= sM->normalizeFFT;
  u_avg[2] /= sM->normalizeFFT;
}


void SpectralAnalysis::_fftw2cub() const
{
  const size_t NlocBlocks = sM->local_infos.size();
  const Real * const data_u = sM->data_u;
  const Real * const data_v = sM->data_v;
  const Real * const data_w = sM->data_w;
  const Real factor = 1 / sM->normalizeFFT;

  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    FluidBlock& b = *(FluidBlock*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(size_t iz=0; iz< (size_t) FluidBlock::sizeZ; ++iz)
    for(size_t iy=0; iy< (size_t) FluidBlock::sizeY; ++iy)
    for(size_t ix=0; ix< (size_t) FluidBlock::sizeX; ++ix) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      b(ix,iy,iz).tmpU += factor * data_u[src_index];
      b(ix,iy,iz).tmpV += factor * data_v[src_index];
      b(ix,iy,iz).tmpW += factor * data_w[src_index];
    }
  }
}

void SpectralAnalysis::run()
{
  _cub2fftw();
  sM->runFwd();
  sM->_compute_analysis();
  sM->runBwd();
  _fftw2cub();

  std::stringstream ssR; ssR<<std::setfill('0')<<std::setw(9)<<sM->sim.step;
  const auto nameV = StreamerVelocityVector::prefix() + ssR.str();
  const auto nameO = StreamerTmpVector::prefix() + ssR.str();
  DumpHDF5_MPI<StreamerVelocityVector, DumpReal>(
    *sM->sim.grid, sM->sim.time, nameV, sM->sim.path4serialization);
  DumpHDF5_MPI<StreamerTmpVector, DumpReal>(
    *sM->sim.grid, sM->sim.time, nameO, sM->sim.path4serialization);
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
    << std::setw(15) << sM->maxGridL
    << " #simulation box length" << "\n";

  f << std::left << std::setw(15) << "tke"
    << std::setw(15) << sM->stats.tke
    << " #turbulent kinetic energy" << "\n";

  f << std::left << std::setw(15) << "eps"
    << std::setw(15) << sM->stats.eps
    << " #Viscous dissipation rate" << "\n";

  f << std::left << std::setw(15) << "eps_f"
    << std::setw(15) << sM->sim.dissipationRate
    << " #Total dissipation rate" << "\n";

  f << std::left << std::setw(15) << "lambda"
    << std::setw(15) << sM->stats.lambda
    << " #Taylor microscale" << "\n";

  f << std::left << std::setw(15) << "Re_lambda"
    << std::setw(15) << sM->stats.Re_lambda
    << " #Turbulent Reynolds based on lambda" << "\n";

  f << std::left << std::setw(15) << "eta"
    << std::setw(15) << std::pow( pow3(sM->sim.nu) / sM->stats.eps, 1.0/4)
    << " #Kolmogorov length scale" << "\n";

  f << std::left << std::setw(15) << "tau_integral"
    << std::setw(15) << sM->stats.tau_integral
    << " #Integral time scale" << "\n";

  f << std::left << std::setw(15) << "tau_eta"
    << std::setw(15) << std::sqrt( sM->sim.nu / sM->stats.eps )
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

  const int nBins = sM->stats.nBin;
  for (int i = 0; i < nBins; i++)
    f << std::left << std::setw(15) << sM->stats.k_msr[i]
                   << std::setw(15) << sM->stats.E_msr[i] << "\n";
  f.flush();
  f.close();

  /*
    std::stringstream ssR_cs2;
    ssR_cs2 << "analysis/spectralAnalysisCs2_" << std::setfill('0')
            << std::setw(9) << nFile;
    f.open(ssR_cs2.str());
    f << std::left << "Cs2 spectrum :" << "\n";
    f << std::left << std::setw(15) << "k * (lBox/2pi)"
                   << std::setw(15) << "Cs2_k" << "\n";
    for (int i = 0; i < nBins; i++)
      f << std::left << std::setw(15) << sM->stats.k_msr[i]
                     << std::setw(15) << sM->stats.cs2_msr[i] << "\n";
    f.flush();
    f.close();
  */
}

void SpectralAnalysis::reset()
{
}

SpectralAnalysis::~SpectralAnalysis()
{
}

CubismUP_3D_NAMESPACE_END
