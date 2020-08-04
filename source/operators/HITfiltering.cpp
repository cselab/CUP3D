//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "HITfiltering.h"

CubismUP_3D_NAMESPACE_BEGIN using namespace cubism;


struct FilteredQuantities
{
  struct Fields
  {
    mutable Real Cs2 = 0;
    mutable Real u  = 0, v  = 0, w  = 0;
    mutable Real uu = 0, uv = 0, uw = 0;
    mutable Real vv = 0, vw = 0, ww = 0;
    mutable Real normalization = 0;

    void reset() const {
      u  = 0; v  = 0; w  = 0; uu = 0; uv = 0; uw = 0; vv = 0; vw = 0; ww = 0;
      normalization = 0;
    }

    void add(const Real U, const Real V, const Real W, const Real fac) const {
      u  += fac * U;   v  += fac * V;   w  += fac * W;
      uu += fac * U*U; uv += fac * U*V; uw += fac * U*W;
      vv += fac * V*V; vw += fac * V*W; ww += fac * W*W;
      normalization += fac;
    }

    void reduce(const Fields & other) const {
      u  += other.u;  v  += other.v;  w  += other.w;
      uu += other.uu; uv += other.uv; uw += other.uw;
      vv += other.vv; vw += other.vw; ww += other.ww;
      normalization += other.normalization;
    }
  };

  const int BPDX, BPDY, BPDZ;
  using data_t = std::vector<std::vector<std::vector<Fields>>>;
  const data_t F = data_t (
      BPDZ, std::vector<std::vector<Fields>>(BPDY, std::vector<Fields>(BPDX)));

  FilteredQuantities(int NX, int NY, int NZ) : BPDX(NX), BPDY(NY), BPDZ(NZ) { }

  const Fields & operator() (const int bx, const int by, const int bz)
  {
    return F[(bz+BPDZ) % BPDZ][(by+BPDY) % BPDY][(bx+BPDX) % BPDX];
  }

  void reset()
  {
    for(int z=0; z<BPDZ; ++z)
      for(int y=0; y<BPDY; ++y)
        for(int x=0; x<BPDX; ++x) F[z][y][x].reset();
  }

  void add(const int bx, const int by, const int bz,
           const Real u, const Real v, const Real w, const Real fac)
  {
    (*this)(bx, by, bz).add(u, v, w, fac);
  }

  void reduce(const FilteredQuantities & other)
  {
    #pragma omp parallel for schedule(static) collapse(3)
    for(int z=0; z<BPDZ; ++z)
      for(int y=0; y<BPDY; ++y)
        for(int x=0; x<BPDX; ++x)
          F[z][y][x].reduce( other.F[z][y][x] );
  }

  void computeCS(const Real h)
  {
    #pragma omp parallel for schedule(static) collapse(3)
    for(int bz=0; bz<BPDZ; ++bz)
    for(int by=0; by<BPDY; ++by)
    for(int bx=0; bx<BPDX; ++bx) {
      (*this)(bx,by,bz).u  /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).v  /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).w  /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).uu /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).uv /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).uw /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).vv /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).vw /= (*this)(bx,by,bz).normalization;
      (*this)(bx,by,bz).ww /= (*this)(bx,by,bz).normalization;
    }

    const Real H = FluidBlock::sizeZ * h;
    #pragma omp parallel for schedule(static) collapse(3)
    for(int bz=0; bz<BPDZ; ++bz)
    for(int by=0; by<BPDY; ++by)
    for(int bx=0; bx<BPDX; ++bx) {
      const Real u  = (*this)(bx, by, bz).u;
      const Real v  = (*this)(bx, by, bz).v;
      const Real w  = (*this)(bx, by, bz).w;
      const Real uu = (*this)(bx, by, bz).uu;
      const Real uv = (*this)(bx, by, bz).uv;
      const Real uw = (*this)(bx, by, bz).uw;
      const Real vv = (*this)(bx, by, bz).vv;
      const Real vw = (*this)(bx, by, bz).vw;
      const Real ww = (*this)(bx, by, bz).ww;
      const Real dudx = (*this)(bx+1,by,  bz  ).u - (*this)(bx-1,by,  bz  ).u;
      const Real dudy = (*this)(bx,  by+1,bz  ).u - (*this)(bx,  by-1,bz  ).u;
      const Real dudz = (*this)(bx,  by,  bz+1).u - (*this)(bx,  by,  bz-1).u;
      const Real dvdx = (*this)(bx+1,by,  bz  ).v - (*this)(bx-1,by,  bz  ).v;
      const Real dvdy = (*this)(bx,  by+1,bz  ).v - (*this)(bx,  by-1,bz  ).v;
      const Real dvdz = (*this)(bx,  by,  bz+1).v - (*this)(bx,  by,  bz-1).v;
      const Real dwdx = (*this)(bx+1,by,  bz  ).w - (*this)(bx-1,by,  bz  ).w;
      const Real dwdy = (*this)(bx,  by+1,bz  ).w - (*this)(bx,  by-1,bz  ).w;
      const Real dwdz = (*this)(bx,  by,  bz+1).w - (*this)(bx,  by,  bz-1).w;

      const Real traceT = (uu - u*u + vv - v*v + ww - w*w) / 3;
      const Real Txx = uu - u*u - traceT;
      const Real Txy = uv - u*v;
      const Real Txz = uw - u*w;
      const Real Tyy = vv - v*v - traceT;
      const Real Tyz = vw - v*w;
      const Real Tzz = ww - w*w - traceT;
      const Real traceS =  (dudx + dvdy + dwdz) / 3;
      const Real Sxx = 0.5/H * (dudx - traceS);
      const Real Sxy = 0.5/H * (dudy + dvdx) / 2;
      const Real Sxz = 0.5/H * (dudz + dwdx) / 2;
      const Real Syy = 0.5/H * (dvdy - traceS);
      const Real Syz = 0.5/H * (dwdy + dvdz) / 2;
      const Real Szz = 0.5/H * (dwdz - traceS);
      Real SS = Sxx*Sxx + Syy*Syy + Szz*Szz + 2*Sxy*Sxy + 2*Sxz*Sxz + 2*Syz*Syz;
      Real TS = Txx*Sxx + Tyy*Syy + Tzz*Szz + 2*Txy*Sxy + 2*Txz*Sxz + 2*Tyz*Syz;
      (*this)(bx, by, bz).Cs2 = - TS / ( 2 * H*H * std::sqrt(2 * SS) * SS );
    }
  }
};

HITfiltering::HITfiltering(SimulationData& s) : Operator(s) {}

void HITfiltering::operator()(const double dt)
{
  if ( not (sim.timeAnalysis>0 && (sim.time+dt) >= sim.nextAnalysisTime) )
    return;

  sim.startProfiler("HITfiltering Kernel");
  const int BPDX = sim.grid->getBlocksPerDimension(0);
  const int BPDY = sim.grid->getBlocksPerDimension(1);
  const int BPDZ = sim.grid->getBlocksPerDimension(2);
  FilteredQuantities filtered(BPDX, BPDY, BPDZ);

  static constexpr int NB = CUP_BLOCK_SIZE;

  #pragma omp parallel for schedule(static) collapse(3)
  for(int biz=0; biz<BPDZ; ++biz)
  for(int biy=0; biy<BPDY; ++biy)
  for(int bix=0; bix<BPDX; ++bix)
  {
    for (int iz = 0; iz < NB; ++iz)
    for (int iy = 0; iy < NB; ++iy)
    for (int ix = 0; ix < NB; ++ix)
    {
      #if 0
        const int bid = bix + biy * BPDX + biz * BPDX * BPDY;
        FluidBlock & block = * (FluidBlock*) vInfo[bid].ptrBlock;
        const Real u = block(ix,iy,iz).u;
        const Real v = block(ix,iy,iz).v;
        const Real w = block(ix,iy,iz).w;
        filtered.add(bix, biy, biz, u, v, w, 1);
      #else
        // linear interp betwen element's block (bix, biy, biz) and second
        // nearest. figure out which from element's index (ix, iy, iz) in block:
        const int nbix = ix >= NB/2 ? bix - 1 : bix + 1;
        const int nbiy = iy >= NB/2 ? biy - 1 : biy + 1;
        const int nbiz = iz >= NB/2 ? biz - 1 : biz + 1;
        // distance from second nearest block along its direction:
        const Real dist_nbix = ix < NB/2 ? NB/2 + ix + 0.5 : 3*NB/2 - ix - 0.5;
        const Real dist_nbiy = iy < NB/2 ? NB/2 + iy + 0.5 : 3*NB/2 - iy - 0.5;
        const Real dist_nbiz = iz < NB/2 ? NB/2 + iz + 0.5 : 3*NB/2 - iz - 0.5;
        // distance from block's center:
        const Real dist_bix = std::fabs(ix + 0.5 - NB/2);
        const Real dist_biy = std::fabs(iy + 0.5 - NB/2);
        const Real dist_biz = std::fabs(iz + 0.5 - NB/2);

        for(int dbz = 0; dbz < 2; ++dbz) // 0 is current block, 1 is
        for(int dby = 0; dby < 2; ++dby) // nearest along z, y, x
        for(int dbx = 0; dbx < 2; ++dbx)
        {
          const int bidx = dbx? nbix : bix;
          const int bidy = dby? nbiy : biy;
          const int bidz = dbz? nbiz : biz;
          const int bid = ( (bidx+BPDX) % BPDX )
                        + ( (bidy+BPDY) % BPDY ) * BPDX
                        + ( (bidz+BPDZ) % BPDZ ) * BPDX * BPDY;
          const BlockInfo & info = vInfo[bid];
          FluidBlock& block = * (FluidBlock*) info.ptrBlock;
          const Real u = block(ix,iy,iz).u;
          const Real v = block(ix,iy,iz).v;
          const Real w = block(ix,iy,iz).w;
          const Real distx = (dbx? dist_nbix : dist_bix) / NB;
          const Real disty = (dby? dist_nbiy : dist_biy) / NB;
          const Real distz = (dbz? dist_nbiz : dist_biz) / NB;
          assert(distx < 1.0 and disty < 1.0 and disty < 1.0);
          //const Real dist =std::sqrt(distx*distx + disty*disty + distz*distz);
          //const Real weight = std::max(1 - dist, (Real) 0);
          const Real weight = (1 - distx) * (1 - disty) * (1 - distz);
          filtered.add(bix, biy, biz, u, v, w, weight);
        }
      #endif
    }
  }

  filtered.computeCS(vInfo[0].h_gridpoint);
  static constexpr Real maxCS2 = 0.15;
  static constexpr Real minCS2 = -0.1;
  static constexpr int nBins = 200;
  int H[nBins] = {0};

  const auto CStoBinID = [&] (const Real CS2) {
    const int signedID = (CS2 - minCS2) * nBins / (maxCS2 - minCS2);
    //printf("%d %e\n", signedID, CS2);
    return std::max((int) 0, std::min(signedID, nBins-1));
  };

  //#pragma omp parallel for reduction(+ : H[:nBins]) schedule(static) collapse(3)
  for(int z=0; z<BPDZ; ++z)
  for(int y=0; y<BPDY; ++y)
  for(int x=0; x<BPDX; ++x) ++ H[ CStoBinID( filtered(x, y, z).Cs2 ) ];

  MPI_Allreduce(MPI_IN_PLACE, H, nBins, MPI_INT, MPI_SUM, sim.app_comm);
  const size_t normalize =  FluidBlock::sizeX * (size_t) sim.bpdx
                          * FluidBlock::sizeY * (size_t) sim.bpdy
                          * FluidBlock::sizeZ * (size_t) sim.bpdz;
  if(sim.rank==0 and not sim.muteAll)
  {
    std::vector<double> buf(nBins, 0);
    for (int i = 0; i < nBins; ++i) buf.push_back( H[i] / (double) normalize);
    FILE * pFile = fopen ("dnsAnalysis.raw", "ab");
    fwrite (buf.data(), sizeof(double), buf.size(), pFile);
    fflush(pFile); fclose(pFile);
  }

  sim.stopProfiler();

  check("SGS");
}

StructureFunctions::StructureFunctions(SimulationData& s) : Operator(s)
{
  std::random_device rd;
  gen.seed(rd());
}

inline Real periodic_distance(const Real x1, const Real x0, const Real extent)
{
  const Real dx = x1 - x0;
  if      (dx >   extent/2) return dx - extent;
  else if (dx <= -extent/2) return dx + extent;
  else return dx;
}
inline Real periodic_distance(const std::array<Real,3> & p1,
                     const std::array<Real,3> & p0,
                     const std::array<Real,3> & extent)
{
  const Real dx = periodic_distance(p1[0], p0[0], extent[0]);
  const Real dy = periodic_distance(p1[1], p0[1], extent[1]);
  const Real dz = periodic_distance(p1[2], p0[2], extent[2]);
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void StructureFunctions::operator()(const double dt)
{
  if ( not (sim.timeAnalysis>0 && (sim.time+dt) >= sim.nextAnalysisTime) )
    return;

  sim.startProfiler("StructureFunctions Kernel");

  std::uniform_int_distribution<size_t> distrib_block(0, vInfo.size()-1);
  std::uniform_int_distribution<int> distrib_elem(0, CUP_BLOCK_SIZE-1);

  const size_t ref_bid = distrib_block(gen);
  const int ref_iz = distrib_elem(gen);
  const int ref_iy = distrib_elem(gen);
  const int ref_ix = distrib_elem(gen);
  const BlockInfo & ref_info = vInfo[ref_bid];
  FluidBlock & ref_block = * (FluidBlock*) ref_info.ptrBlock;
  const FluidElement & ref_elem = ref_block(ref_ix, ref_iy, ref_iz);
  const std::array<Real,3> ref_pos = ref_info.pos<Real>(ref_ix, ref_iy, ref_iz);

  const Real h, delta_h = sim.uniformH(), 2 * h;
  const size_t n_shells = sim.maxextent / delta_h;

  std::vector<double> sum_S2(n_shells, 0), sumsq_S2(n_shells, 0);
  std::vector<double> sum_S3(n_shells, 0), sumsq_S3(n_shells, 0);
  std::vector<int>  counts(n_shells, 0);

  #pragma omp parallel for schedule(static) reduction(+ : counts[:n_shells],   \
                                                          sum_S2[:n_shells],   \
                                                          sumsq_S2[:n_shells], \
                                                          sum_S3[:n_shells],   \
                                                          sumsq_S3[:n_shells])
  for (size_t i = 0; i < vInfo.size(); ++i)
  {
    const BlockInfo & info = vInfo[i];
    FluidBlock& block = * (FluidBlock*) info.ptrBlock;

    for (int iz = 0; iz < CUP_BLOCK_SIZE; ++iz)
    for (int iy = 0; iy < CUP_BLOCK_SIZE; ++iy)
    for (int ix = 0; ix < CUP_BLOCK_SIZE; ++ix)
    {
      const std::array<Real,3> pos = info.pos<Real>(ix, iy, iz);
      const Real delta = periodic_distance(pos, ref_pos, sim.extent);
      const size_t shell_id = delta / delta_h;
      if (shell_id >= n_shells) continue;
      const Real S2_ui = std::pow(block(ix,iy,iz).u - ref_elem.u, 2);
      const Real S2_vi = std::pow(block(ix,iy,iz).v - ref_elem.v, 2);
      const Real S2_wi = std::pow(block(ix,iy,iz).w - ref_elem.w, 2);
      const Real S3_ui = std::pow(block(ix,iy,iz).u - ref_elem.u, 3);
      const Real S3_vi = std::pow(block(ix,iy,iz).v - ref_elem.v, 3);
      const Real S3_wi = std::pow(block(ix,iy,iz).w - ref_elem.w, 3);
      sum_S2  [shell_id] += S2_ui + S2_vi + S2_wi;
      sumsq_S2[shell_id] += S2_ui*S2_ui + S2_vi*S2_vi + S2_wi*S2_wi;
      sum_S3  [shell_id] += S3_ui + S3_vi + S3_wi;
      sumsq_S3[shell_id] += S3_ui*S3_ui + S3_vi*S3_vi + S3_wi*S3_wi;
      counts  [shell_id] += 1;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,   sum_S2.data(), n_shells, MPI_DOUBLE, MPI_SUM, sim.app_comm);
  MPI_Allreduce(MPI_IN_PLACE, sumsq_S2.data(), n_shells, MPI_DOUBLE, MPI_SUM, sim.app_comm);
  MPI_Allreduce(MPI_IN_PLACE,   sum_S3.data(), n_shells, MPI_DOUBLE, MPI_SUM, sim.app_comm);
  MPI_Allreduce(MPI_IN_PLACE, sumsq_S3.data(), n_shells, MPI_DOUBLE, MPI_SUM, sim.app_comm);
  MPI_Allreduce(MPI_IN_PLACE,   counts.data(), n_shells, MPI_INT,    MPI_SUM, sim.app_comm);

  if(sim.rank==0 and not sim.muteAll)
  {
    std::vector<double> buffer = sum_S2;
    buffer.insert(buffer.end(), sumsq_S2.begin(), sumsq_S2.end());
    buffer.insert(buffer.end(),   sum_S3.begin(),   sum_S3.end());
    buffer.insert(buffer.end(), sumsq_S3.begin(), sumsq_S3.end());
    buffer.insert(buffer.end(),   counts.begin(),   counts.end());
    FILE * pFile = fopen ("structureFunctionsAnalysis.raw", "ab");
    fwrite (buf.data(), sizeof(double), buf.size(), pFile);
    fflush(pFile); fclose(pFile);
  }

  sim.stopProfiler();

  check("SGS");
}

CubismUP_3D_NAMESPACE_END
