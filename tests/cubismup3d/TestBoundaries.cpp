#include "Utils.h"
#include "../../source/Simulation.h"
#include "../../source/operators/CellwiseOperator.h"

using namespace cubism;
using namespace cubismup3d;

static constexpr int CELLS_X = 16*8;
static constexpr int CELLS_Y = 4*8;
static constexpr int CELLS_Z = 4*8;

/* Expected value of a given cell. */
static double getElementValue(int abs_ix, int abs_iy, int abs_iz)
{
  return 0.1 + abs_ix + 1000 * abs_iy + 1000000 * abs_iz;
}

/* Test one stencil on periodic boundaries with a ready simulation object. */
void _testPeriodicBoundaries(Simulation &s, StencilInfo stencil, int dx, int dy, int dz)
{
  // Reset the grid to some initial vlaue.
  applyKernel(s.sim, [](CellInfo info, FluidElement &e) {
    e.u = getElementValue(info.get_abs_ix(), info.get_abs_iy(), info.get_abs_iz());
  });

  // Shift everything in the direction (dx, dy, dz);
  applyStencilKernel(
    s.sim,
    stencil,
    [dx, dy, dz](const StencilKernelLab &lab, CellInfo info, FluidElement &out) {
      out.tmpU = lab(-dx, -dy, -dz).u;
    }
  );

  // Check if everything is correct now.
  applyKernel(s.sim, [dx, dy, dz](CellInfo info, FluidElement &e) {
    int aix = (info.get_abs_ix() - dx + CELLS_X) % CELLS_X;
    int aiy = (info.get_abs_iy() - dy + CELLS_Y) % CELLS_Y;
    int aiz = (info.get_abs_iz() - dz + CELLS_Z) % CELLS_Z;
    double expected = getElementValue(aix, aiy, aiz);
    if (expected != e.tmpU) {
      fprintf(stderr, "local=(%d %d %d) global=(%d %d %d) block=(%d %d %d) is %f instead of %f\n",
              info.ix, info.iy, info.iz,
              info.get_abs_ix(),
              info.get_abs_iy(),
              info.get_abs_iz(),
              info.block_info.index[0],
              info.block_info.index[1],
              info.block_info.index[2],
              e.tmpU, expected);
      fprintf(stderr, "Failed at shift (%d, %d, %d)\n", dx, dy, dz);;
      exit(1);
    }
  });
}

/* Test stencils on periodic boundaries. */
bool testPeriodicBoundaries()
{
  const std::vector<std::string> argv{
    "-cfl", "0.1",
    "-BC_x", "periodic",
    "-BC_y", "periodic",
    "-BC_z", "periodic",
    "-bpdx", computeNumBlocksArg(CELLS_X, FluidBlock::sizeX),
    "-bpdy", computeNumBlocksArg(CELLS_Y, FluidBlock::sizeY),
    "-bpdz", computeNumBlocksArg(CELLS_Z, FluidBlock::sizeZ),
  };
  auto s = createSimulation(MPI_COMM_WORLD, argv);

  // Try out 3 different stencils.
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), -1, 0, 0);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), +1, 0, 0);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), 0, -1, 0);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), 0, +1, 0);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), 0, 0, -1);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), 0, 0, +1);
  _testPeriodicBoundaries(*s, StencilInfo(-1, -1, -1, 2, 2, 2, false, {{FE_U}}), +1, 0, 0);
  _testPeriodicBoundaries(*s, StencilInfo(-2, -2, -2, 3, 3, 3, true,  {{FE_U}}), -1, +1, +2);

  return true;
}

int main(int argc, char **argv)
{
  initMPI(&argc, &argv);

  CUP_RUN_TEST(testPeriodicBoundaries);

  finalizeMPI();
}
