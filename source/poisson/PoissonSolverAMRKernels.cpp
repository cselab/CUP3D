#include "PoissonSolverAMRKernels.h"

#include <cassert>

/*
Optimization comments:
  - The innermost loop has to be very simple in order for the compiler to
    optimize it. Temporary accumulators and storage arrays have to be used to
    enable vectorization.

  - In order to vectorize the stencil, the shifted west and east pointers to p
    have to be provided separately, without the compiler knowing that they are
    related to the same buffer.

  - The same would be true for south, north, back and front shifts, but we pad
    the p block not with +/-1 padding but with +/-4, to have everything nicely
    aligned with respect to the 32B boundary. This was tested only on AVX-256.
    However, the intentional offset of 32B should be favorable for AVX-512 as
    well.

  - For correctness, the p pointers must not have __restrict__, since p,
    pW and pE do overlap. (Not important here though, since we removed the
    convergence for loop, see below). All other arrays do have __restrict__, so
    this does not affect vectorization anyway.

  - The outer for loop that repeats the kernel until convergence breaks the
    vectorization of the stencil in gcc, hence it was removed from this file.

  - Putting this loop into another function does not help, since the compiler
    merges the two functions and breaks the vectorization. This can be fixed by
    adding `static __attribute__((noinline))` to the kernel function, but it is
    a bit risky, and doesn't seem to improve the generated code. The cost of
    the bare function call here is about 3ns.

  - Not tested here, but unaligned access can be up to 2x slower than aligned,
    so it is important to ensure alignment.
    https://www.agner.org/optimize/blog/read.php?i=423


Compilation hints:
  - If gcc is used, the Ax--p stencil won't be vectorized unless version 11 or
    later is used.

  - -ffast-math might affect ILP and reductions. Not verified.

  - Changing the order of operations may cause the compiler to produce
    different operation order and hence cause the number of convergence
    iterations to change.

  - To show the assembly, use e.g.
      objdump -dS -Mintel --no-show-raw-insn PoissonSolverAMRKernels.cpp.o > PoissonSolverAMRKErnels.cpp.lst

  - With gcc 11, it might be necessary to use "-g -gdwarf-4" instead of "-g"
    for objdump to work. For more information look here:
    https://gcc.gnu.org/gcc-11/changes.html


Benchmarks info for Broadwell CPU with AVX2 (256-bit):
  - Computational limit is 2x 256-bit SIMD FMAs per cycle == 16 FLOPs/cycle.

  - Memory limit for L1 cache is 2x256-bit reads and 1x256-bit write per cycle.
    See "Haswell and Broadwell pipeline", section "Read and write bandwidth":
    https://www.agner.org/optimize/microarchitecture.pdf

    These amount to 64B reads and 32B writes per cycle, however we get about
    80% of that, consistent with benchmarks here:
    https://www.agner.org/optimize/blog/read.php?i=423

  - The kernels below are memory bound.
*/

namespace cubismup3d {
namespace poisson_kernels {

/// Update `r -= a * Ax` and return `sum(r^2)`.
static inline Real subAndSumSqr(
    Block & __restrict__ r_,
    const Block & __restrict__ Ax_,
    Real a)
{
  // Equivalent 3D code:
  // Real s = {};
  // for (int iz = 0; iz < NZ; ++iz)
  // for (int iy = 0; iy < NY; ++iy)
  // for (int ix = 0; ix < NX; ++ix) {
  //     r_[iz][iy][ix] -= a * Ax_[iz][iy][ix];
  //     s += r_[iz][iy][ix] * r_[iz][iy][ix];
  // }
  // return s;

  // The block structure is not important here, we can treat it as a contiguous
  // array. However, we group into groups of length 16, to help with ILP and
  // vectorization.
  constexpr int MX = 16;
  constexpr int MY = NX * NY * NZ / MX;
  using SquashedBlock = Real[MY][MX];
  static_assert(NX * NY % 16 == 0 && sizeof(Block) == sizeof(SquashedBlock));
  SquashedBlock & __restrict__ r = (SquashedBlock &)r_;
  SquashedBlock & __restrict__ Ax = (SquashedBlock &)Ax_;

  // This kernel does not reach neither the compute nor the memory bound.
  // The problem could be high latency of FMA instructions.
  Real s[MX] = {};
  for (int jy = 0; jy < MY; ++jy) {
    for (int jx = 0; jx < MX; ++jx)
      r[jy][jx] -= a * Ax[jy][jx];
    for (int jx = 0; jx < MX; ++jx)
      s[jx] += r[jy][jx] * r[jy][jx];
  }
  return sum(s);
}

template <typename T>
static inline T *assumeAligned(T *ptr, unsigned align, unsigned offset = 0)
{
  if (sizeof(Real) == 8 || sizeof(Real) == 4) {
    // if ((uintptr_t)ptr % align != offset)
    //   throw std::runtime_error("alignment failed");
    assert((uintptr_t)ptr % align == offset);

    // Works with gcc, clang and icc.
    return (T *)__builtin_assume_aligned(ptr, align, offset);
  } else {
    return ptr;  // No alignment for long long.
  }
};

Real kernelPoissonGetZInner(
    PaddedBlock &p_,
    const Real *pW_,
    const Real *pE_,
    Block & __restrict__ Ax_,
    Block & __restrict__ r_,
    Block & __restrict__ block_,
    const Real sqrNorm0,
    const Real rr)
{
  PaddedBlock &p = *assumeAligned(&p_, 64, 64 - xShift * sizeof(Real));
  const PaddedBlock &pW = *(PaddedBlock *)pW_;  // Aligned to 64B + 24 (for doubles).
  const PaddedBlock &pE = *(PaddedBlock *)pE_;  // Aligned to 64B + 40 (for doubles).
  Block & __restrict__ Ax = *assumeAligned(&Ax_, 64);
  Block & __restrict__ r = *assumeAligned(&r_, 64);
  Block & __restrict__ block = *assumeAligned(&block_, kBlockAlignment);

  // Broadwell: 6.0-6.6 FLOP/cycle, depending probably on array alignments.
  Real a2Partial[NX] = {};
  for (int iz = 0; iz < NZ; ++iz)
  for (int iy = 0; iy < NY; ++iy) {
    // On Broadwell and earlier it might be beneficial to turn some of these
    // a+b additions into FMAs of form 1*a+b, because those CPUs can do 2
    // FMAs/cycle and only 1 ADD/cycle. However, it wouldn't be simple to
    // convience the compiler to do so, and it wouldn't matter from Skylake on.
    // https://www.agner.org/optimize/blog/read.php?i=415

    Real tmpAx[NX];
    for (int ix = 0; ix < NX; ++ix) {
      tmpAx[ix] = pW[iz + 1][iy + 1][ix + xShift]
                + pE[iz + 1][iy + 1][ix + xShift]
                - 6.0 * p[iz + 1][iy + 1][ix + xShift];
    }

    // This kernel is memory bound. After unrolling the iy loop the compiler
    // should figure out that some loads can be reused between consecutive iy.

    // You might be tempted to merge these two loops (i.e. to preserve
    // symmetry), but that kills the vectorization in gcc 11. At least only
    // when the `p = r + beta * p` loop is enabled...
    for (int ix = 0; ix < NX; ++ix)
      tmpAx[ix] += p[iz + 1][iy][ix + xShift];
    for (int ix = 0; ix < NX; ++ix)
      tmpAx[ix] += p[iz + 1][iy + 2][ix + xShift];

    for (int ix = 0; ix < NX; ++ix)
      tmpAx[ix] += p[iz][iy + 1][ix + xShift];
    for (int ix = 0; ix < NX; ++ix)
      tmpAx[ix] += p[iz + 2][iy + 1][ix + xShift];

    for (int ix = 0; ix < NX; ++ix)
      Ax[iz][iy][ix] = tmpAx[ix];

    for (int ix = 0; ix < NX; ++ix)
      a2Partial[ix] += p[iz + 1][iy + 1][ix + xShift] * tmpAx[ix];
  }
  const Real a2 = sum(a2Partial);
  const Real a = rr / (a2 + 1e-55);

  // In theory, merging this memory bound kernel with the next one should be
  // beneficial since the other one is not memory bound. But, merging seems to
  // degrade the performance.
  for (int iz = 0; iz < NZ; ++iz)
  for (int iy = 0; iy < NY; ++iy)
  for (int ix = 0; ix < NX; ++ix)
    block[iz][iy][ix] += a * p[iz + 1][iy + 1][ix + xShift];

  // r -= a * Ax, sqrSum = sum(r^2)
  const Real sqrSum = subAndSumSqr(r, Ax, a);

  const Real beta = sqrSum / (rr + 1e-55);
  const Real sqrNorm = (Real)1 / (N * N) * sqrSum;

  if (sqrNorm < 1e-14 * sqrNorm0 || sqrNorm < 1e-32)
    return (Real)0.0;

  // Kernel: 2 reads + 1 write + 2 FLOPs per cell -> limit is L1 cache.
  // Broadwell: 6.5 FLOP/cycle, 52+26 B/cycle
  for (int iz = 0; iz < NZ; ++iz)
  for (int iy = 0; iy < NY; ++iy)
  for (int ix = 0; ix < NX; ++ix) {
    p[iz + 1][iy + 1][ix + xShift] =
        r[iz][iy][ix] + beta * p[iz + 1][iy + 1][ix + xShift];
  }

  const Real rrNew = sqrSum;
  return rrNew;
}

}  // namespace poisson_kernels
}  // namespace cubismup3d
