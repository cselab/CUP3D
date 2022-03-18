#include "../../source/Simulation.h"
#include "CellwiseOperator.h"
#include "Utils.h"

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
static void _testPeriodicBoundaries(Simulation &s, StencilInfo stencil, int dx, int dy, int dz)
{
  // Reset the grid to some initial vlaue.
  applyKernel(s.sim, [](CellInfo info, FluidElement &e) {
    e.u = getElementValue(info.getAbsIX(), info.getAbsIY(), info.getAbsIZ());
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
    int aix = (info.getAbsIX() - dx + CELLS_X) % CELLS_X;
    int aiy = (info.getAbsIY() - dy + CELLS_Y) % CELLS_Y;
    int aiz = (info.getAbsIZ() - dz + CELLS_Z) % CELLS_Z;
    double expected = getElementValue(aix, aiy, aiz);
    if (expected != e.tmpU) {
      fprintf(stderr, "local=(%d %d %d) global=(%d %d %d) block=(%d %d %d) is %f instead of %f\n",
              info.ix, info.iy, info.iz,
              info.getAbsIX(),
              info.getAbsIY(),
              info.getAbsIZ(),
              info.blockInfo.index[0],
              info.blockInfo.index[1],
              info.blockInfo.index[2],
              e.tmpU, expected);
      fprintf(stderr, "Failed at shift (%d, %d, %d)\n", dx, dy, dz);;
      exit(1);
    }
  });
}

/* Test stencils on periodic boundaries. */
static void testPeriodicBoundaries(int levelMax, int levelStart)
{
  const std::vector<std::string> argv{
    "-cfl", "0.1",
    "-BC_x", "periodic",
    "-BC_y", "periodic",
    "-BC_z", "periodic",
    "-bpdx", computeNumBlocksArg(CELLS_X, FluidBlock::sizeX),
    "-bpdy", computeNumBlocksArg(CELLS_Y, FluidBlock::sizeY),
    "-bpdz", computeNumBlocksArg(CELLS_Z, FluidBlock::sizeZ),
    "-levelMax", std::to_string(levelMax),
    "-levelStart", std::to_string(levelStart),
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
}

int main(int argc, char **argv)
{
  int levelMax = 1;
  // TODO: Test for levelMax >= 2 (requires AMR-aware CellInfo).
  /*
  if (argc >= 3) {
    fprintf(stderr, "usage: %s [levelMax]\n", argv[0]);
    return 1;
  }
  if (argc >= 2 && 1 != sscanf(argv[1], "%d", &levelMax)) {
    fprintf(stderr, "expected int, got '%s'\n", argv[1]);
    return 1;
  }
  */

  initMPI(&argc, &argv);

  testPeriodicBoundaries(levelMax, levelMax - 1);

  finalizeMPI();
}
