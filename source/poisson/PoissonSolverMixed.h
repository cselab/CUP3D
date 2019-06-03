//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PoissonSolverMixed_h
#define CubismUP_3D_PoissonSolverMixed_h

#include "poisson/PoissonSolver.h"

CubismUP_3D_NAMESPACE_BEGIN

class PoissonSolverMixed : public PoissonSolver
{
  void * fwd, * bwd;
  ptrdiff_t alloc_local=0,local_n0=0,local_0_start=0,local_n1=0,local_1_start=0;
  const double h = sim.uniformH();
  inline bool DFT_X() const { return sim.BCx_flag == periodic; }
  inline bool DFT_Y() const { return sim.BCy_flag == periodic; }
  inline bool DFT_Z() const { return sim.BCz_flag == periodic; }

 protected:

  template<bool DFTX, bool DFTY, bool DFTZ> void _solve()
  {
    // if BC flag == 1 fourier, else cosine transform
    const Real normX = (DFTX ? 1.0 : 0.5) / ( gsize[0]*h );
    const Real normY = (DFTY ? 1.0 : 0.5) / ( gsize[1]*h );
    const Real normZ = (DFTZ ? 1.0 : 0.5) / ( gsize[2]*h );
    const Real waveFactX = (DFTX ? 2 : 1) * M_PI / ( gsize[0]*h );
    const Real waveFactY = (DFTY ? 2 : 1) * M_PI / ( gsize[1]*h );
    const Real waveFactZ = (DFTZ ? 2 : 1) * M_PI / ( gsize[2]*h );
    const Real norm_factor = normX * normY * normZ;
    Real *const in_out = data;
    const long nKx = static_cast<long>(gsize[0]);
    const long nKy = static_cast<long>(gsize[1]);
    const long nKz = static_cast<long>(gsize[2]);
    const long shifty = static_cast<long>(local_1_start);
    #pragma omp parallel for schedule(static)
    for(long j = 0; j<static_cast<long>(local_n1); ++j)
    for(long i = 0; i<static_cast<long>(gsize[0]); ++i)
    for(long k = 0; k<static_cast<long>(gsize[2]); ++k)
    {
      const size_t linidx = (j*gsize[0] +i)*gsize[2] + k;
      const long J = shifty + j; //memory index plus shift due to decomp
      const long kx = DFTX ? ((i <= nKx/2) ? i : nKx-i) : i;
      const long ky = DFTY ? ((J <= nKy/2) ? J : nKy-J) : J;
      const long kz = DFTZ ? ((k <= nKz/2) ? k : nKz-k) : k;
      const Real rkx = ( kx + (DFTX ? 0 : (Real)0.5 ) ) * waveFactX;
      const Real rky = ( ky + (DFTY ? 0 : (Real)0.5 ) ) * waveFactY;
      const Real rkz = ( kz + (DFTZ ? 0 : (Real)0.5 ) ) * waveFactZ;
      in_out[linidx] *= - norm_factor/(rkx*rkx + rky*rky + rkz*rkz);
    }
    //if (shifty==0 && DFTX && DFTY && DFTZ) in_out[0] = 0;
    if (shifty==0) in_out[0] = 0;
  }

 public:

  PoissonSolverMixed(SimulationData & s);

  void solve();

  ~PoissonSolverMixed();
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_PoissonSolverMixed_h
