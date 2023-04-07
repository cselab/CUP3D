#pragma once

#include "../Definitions.h"

namespace cubismup3d {
namespace diffusion_kernels {

static constexpr int NX = ScalarBlock::sizeX;
static constexpr int NY = ScalarBlock::sizeY;
static constexpr int NZ = ScalarBlock::sizeZ;
static constexpr int N = NX * NY * NZ;
using Block = Real[NZ][NY][NX];

// Assuming NX == 8, pad not with +/-1 but with +/-4 to reduce the number of
// unaligned memory accesses, which are up to 2x slower than aligned ones.
static constexpr int xPad = 4;
using PaddedBlock = Real[NZ + 2][NY + 2][NX + 2 * xPad];

template <int N>
static inline Real sum(const Real (&a)[N])
{
  Real s = 0;
  for (int ix = 0; ix < N; ++ix)
    s += a[ix];
  return s;
}

// Simple implementation of the kernel.
Real kernelDiffusionGetZInnerReference(
    PaddedBlock & __restrict__ p_,
    Block & __restrict__ Ax_,
    Block & __restrict__ r_,
    Block & __restrict__ block_,
    const Real sqrNorm0,
    const Real rr);

// Optimized implementation. See the alignment requirements in the code!
Real kernelDiffusionGetZInner(
    PaddedBlock &p,
    const Real *pW,
    const Real *pE,
    Block & __restrict__ Ax,
    Block & __restrict__ r,
    Block & __restrict__ block,
    Real sqrNorm0,
    Real rr, const Real coefficient);

void getZImplParallel(const std::vector<cubism::BlockInfo>& vInfo, const Real nu, const Real dt);

}  // namespace diffusion_kernels
}  // namespace cubismup3d
