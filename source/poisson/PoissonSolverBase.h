#pragma once

#include <memory>
#include "../Definitions.h"

namespace cubismup3d {

struct SimulationData;

class PoissonSolverBase
{
public:
  ~PoissonSolverBase() = default;
  virtual void solve() = 0;
};

std::shared_ptr<PoissonSolverBase> makePoissonSolver(SimulationData& s);
} // namespace cubismup3d
