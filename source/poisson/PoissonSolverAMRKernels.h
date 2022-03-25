#pragma once

#include "../Definitions.h"

namespace cubismup3d {
namespace poisson_kernels {

static constexpr int NX = ScalarBlock::sizeX;
static constexpr int NY = ScalarBlock::sizeY;
static constexpr int NZ = ScalarBlock::sizeZ;
static constexpr int N = NX * NY * NZ;
static constexpr int xShift = 4;
using Block = Real[NZ][NY][NX];

// Pad not with +/-1 but with +/-4 to reduce the number of unaligned memory
// accesses, which are up to 2x slower than aligned ones.
using PaddedBlock = Real[NZ + 2][NY + 2][NX + 8];

template <int N>
static inline Real sum(const Real (&a)[N])
{
  Real s = 0;
  for (int ix = 0; ix < N; ++ix)
    s += a[ix];
  return s;
}

// See the alignment requirements in the code!
Real kernelPoissonGetZInner(
    PaddedBlock &p,
    const Real *pW,
    const Real *pE,
    Block & __restrict__ Ax,
    Block & __restrict__ r,
    Block & __restrict__ block,
    Real sqrNorm0,
    Real rr);

}  // namespace poisson_kernels
}  // namespace cubismup3d
