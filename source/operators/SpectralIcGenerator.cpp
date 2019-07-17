//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "operators/SpectralIcGenerator.h"

#include<random>
#include<iomanip>
#include<sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;


SpectralIcGenerator::SpectralIcGenerator(SimulationData &s)
{
  if (s.spectralManip == nullptr)
    s.spectralManip = new SpectralManip(s);

  sM = s.spectralManip;
  sM->prepareBwd();
}

energySpectrum SpectralIcGenerator::_generateTarget()
{
  std::vector<Real> k; std::vector<Real> E;

  if (sM->sim.spectralIC=="cbc"){
    // Set-up Energy Spectrum from COMTE-BELLOT-CORRSIN experiment :
    // Wave numbers are given from experiment units [cm-1]
    //     => Need to convert to [m^-1] and non-dimensionalize with L_ref
    // Energy is given from experiment units [cm3/s2]
    //     => Need to convert to [m^3/s^2] and non-dimensionalize with L_ref and U0
    k = std::vector<Real> (19, 0.0); E = std::vector<Real> (19, 0.0);
    k[0]  =  0.20;  E[0]  = 129.0;
    k[1]  =  0.25;  E[1]  = 230.0;
    k[2]  =  0.30;  E[2]  = 322.0;
    k[3]  =  0.40;  E[3]  = 435.0;
    k[4]  =  0.50;  E[4]  = 457.0;
    k[5]  =  0.70;  E[5]  = 380.0;
    k[6]  =  1.00;  E[6]  = 270.0;
    k[7]  =  1.50;  E[7]  = 168.0;
    k[8]  =  2.00;  E[8]  = 120.0;
    k[9]  =  2.50;  E[9]  =  89.0;
    k[10] =  3.00;  E[10] =  70.3;
    k[11] =  4.00;  E[11] =  47.0;
    k[12] =  6.00;  E[12] =  24.7;
    k[13] =  8.00;  E[13] =  12.6;
    k[14] = 10.00;  E[14] =  7.42;
    k[15] = 12.50;  E[15] =  3.96;
    k[16] = 15.00;  E[16] =  2.33;
    k[17] = 17.50;  E[17] =  1.34;
    k[18] = 20.00;  E[18] =   0.8;
    const Real M = 5.08;

    const Real u_ref = sqrt(3.0/2)*22.2;
    const Real L_ref = 10.8*M;

    for (std::size_t idx = 0; idx < k.size(); idx++){
      k[idx] *= L_ref * 0.01 * 2 * M_PI;
      E[idx] /= u_ref * u_ref * L_ref;
    }
  }
  else if (sM->sim.spectralIC=="art") {
    // Initial spectrum from :
    // P. Sullivan, Neal & Mahalingam, Shankar & Kerr, Robert. (1994).
    // Deterministic forcing of homogeneous, isotropic turbulence.
    // Physics of Fluids. 6. 10.1063/1.868274.
    const Real k0   = sM->sim.k0;
    const Real tke0 = sM->sim.tke0;
    k = std::vector<Real> (sM->nBin, 0.0); E = std::vector<Real> (sM->nBin, 0.0);
    for (int i=0; i<sM->nBin; i++){
      k[i] = (i+0.5) * sM->binSize;
      E[i] = 32 * tke0 / (3*k0) * sqrt(2.0/M_PI) * pow(k[i]/k0,4.0) * exp(-2*pow(k[i]/k0,2.0));
    }
  }

  else if (sM->sim.spectralIC=="fromFile") {
    std::ifstream inFile;
    std::string fileName = sM->sim.spectralICFile;
    inFile.open(fileName);
    if (!inFile){
      std::cout<<"SpectralICGenerator: cannot open file :"<<fileName<<std::endl;
      abort();
    }
    for(std::string line; std::getline(inFile, line); )
    {
      std::istringstream in(line);
      Real k_r, E_r;
      in >> k_r >> E_r;
      k.push_back(k_r);
      E.push_back(E_r);
    }
    energySpectrum target = energySpectrum(k,E);

    // Set target tke
    Real k_eval = 0.0, tke0 = 0.0;
    for (int i = 0; i<sM->maxGridSize; i++){
      k_eval = (i+1) * 2*M_PI / sM->maxBoxLength;
      tke0  += target.interpE(k_eval);
    }
    sM->sim.tkeTgt = tke0;
  }

  energySpectrum target = energySpectrum(k,E);

  return target;
}

void SpectralIcGenerator::_compute()
{
  std::random_device seed;
  std::mt19937 gen(seed());
  std::uniform_real_distribution<Real> randUniform(0,1);

  energySpectrum target = _generateTarget();

  Real tke = 0.0;
  fft_c *const cplxData_u = (fft_c *) sM->data_u;
  fft_c *const cplxData_v = (fft_c *) sM->data_v;
  fft_c *const cplxData_w = (fft_c *) sM->data_w;
  #pragma omp parallel for reduction(+:tke)
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
    const Real k_xy = sqrt(kx*kx + ky*ky);
    const Real k_norm = sqrt(k2);

    const Real E_k = target.interpE(k_norm);
    const Real amp = (k2==0)? 0. : sqrt(E_k /(2*M_PI* pow(k_norm/sM->waveFactor[0], 2)));

    const Real theta1 = randUniform(gen)*2.0*M_PI;
    const Real theta2 = randUniform(gen)*2.0*M_PI;
    const Real phi    = randUniform(gen)*2.0*M_PI;

    const Real noise_a[2] = {amp * cos(theta1) * cos(phi), amp * sin(theta1) * cos(phi)};
    const Real noise_b[2] = {amp * cos(theta2) * sin(phi), amp * sin(theta2) * sin(phi)};

    const Real fac = k_norm*k_xy;
    const Real invFac = (fac==0)? 0. : 1.0/fac;

    const int mult = (k==0) or (k==sM->nKz/2) ? 1 : 2;

    cplxData_u[linidx][0] = (k_norm==0)? 0.0 : invFac * (noise_a[0] * k_norm * ky + noise_b[0] * kx * kz );
    cplxData_u[linidx][1] = (k_norm==0)? 0.0 : invFac * (noise_a[0] * k_norm * ky + noise_b[1] * kx * kz );

    cplxData_v[linidx][0] = (k_norm==0)? 0.0 : invFac * (noise_b[0] * ky * kz - noise_a[0] * k_norm * kx );
    cplxData_v[linidx][1] = (k_norm==0)? 0.0 : invFac * (noise_b[1] * ky * kz - noise_a[1] * k_norm * kx );

    cplxData_w[linidx][0] = (k_norm==0)? 0.0 : -noise_b[0] * k_xy / k_norm;
    cplxData_w[linidx][1] = (k_norm==0)? 0.0 : -noise_b[1] * k_xy / k_norm;

    tke += 0.5*mult*(pow2_cplx(cplxData_u[linidx]) + pow2_cplx(cplxData_v[linidx]) + pow2_cplx(cplxData_w[linidx]));
  }
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, sM->m_comm);

  if (sM->sim.spectralForcing and sM->sim.tkeTgt==0) sM->sim.tkeTgt = tke;
}

void SpectralIcGenerator::_fftw2cub() const
{
  const size_t NlocBlocks = sM->local_infos.size();
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) sM->local_infos[i].ptrBlock;
    const size_t offset = sM->_offset( sM->local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = sM->_dest(offset, iz, iy, ix);
      b(ix,iy,iz).u = sM->data_u[src_index];
      b(ix,iy,iz).v = sM->data_v[src_index];
      b(ix,iy,iz).w = sM->data_w[src_index];
    }
  }
}

void SpectralIcGenerator::run()
{
  _compute();
  sM->runBwd();
  _fftw2cub();
}

CubismUP_3D_NAMESPACE_END
#undef MPIREAL
