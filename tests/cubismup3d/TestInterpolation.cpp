#include "../Utils.h"
#include "Simulation.h"
#include "operators/CellwiseOperator.h"
#include "operators/LinearInterpolation.h"

double getValue(std::array<double, 3> p) {
  // Must be linear (we are testing linear interpolation).
  return 100.0 * p[0] + 2343 * p[1] + 123. * p[2];
}

bool testLinearInterpolation() {
  constexpr int NUM_POINTS = 10000;
  constexpr double extent = 100.0;

  // Prepare simulation data and the simulation object.
  auto prepareSimulationData = []() {
    SimulationData SD{MPI_COMM_WORLD};
    SD.CFL = 0.1;
    SD.BCx_flag = periodic;  // <--- Periodic boundaries.
    SD.BCy_flag = periodic;
    SD.BCz_flag = periodic;
    SD.extent[0] = extent;
    SD.setCells(64, 64, 64);
    return SD;
  };
  Simulation S{prepareSimulationData()};

  // Reset the grid to some initial vlaue.
  cubismup3d::applyKernel(S.sim, [](cubismup3d::CellInfo info, FluidElement &e) {
    e.u = getValue(info.get_pos());
  });

  // Generate random points. Avoid boundaries.
  std::vector<std::array<double, 3>> points;
  points.reserve(NUM_POINTS);
  for (int i = 0; i < NUM_POINTS; ++i) {
    double x = extent * (0.1 + 0.8 / RAND_MAX * rand());
    double y = extent * (0.1 + 0.8 / RAND_MAX * rand());
    double z = extent * (0.1 + 0.8 / RAND_MAX * rand());
    points.push_back({x, y, z});
  }
  std::vector<double> result(NUM_POINTS);

  // Interpolate.
  cubismup3d::linearCellCenteredInterpolation(
      S.sim,
      points,
      [](const FluidElement &e) { return e.u; },      // What to interpolate.
      [&result](int k, double v) { result[k] = v; },  // Where to store the result.
      std::vector<int>({CUP_ELEMENT_INDEX(u)}));      // Components for stencil.

  // Check if result correct.
  for (int i = 0; i < NUM_POINTS; ++i) {
    double expected = getValue(points[i]);
    if (std::fabs(expected - result[i]) > 1e-9) {
      fprintf(stderr, "Expected %lf, got %lf. Point is (%lf, %lf %lf).\n",
              expected, result[i], points[i][0], points[i][1], points[i][2]);
    }
  }

  return true;
}

int main(int argc, char **argv) {
  cubismup3d::tests::init(&argc, &argv);

  CUP_RUN_TEST(testLinearInterpolation);

  cubismup3d::tests::finalize();
}
