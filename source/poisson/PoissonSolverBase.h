#pragma once

#include <memory>
#include "../Definitions.h"

namespace cubismup3d {

struct SimulationData;

class PoissonSolverBase
{
public:
  virtual ~PoissonSolverBase() = default;
  virtual void solve() = 0;
protected:
  typedef typename ScalarGrid::BlockType BlockType;
};

std::shared_ptr<PoissonSolverBase> makePoissonSolver(SimulationData& s);
} // namespace cubismup3d
