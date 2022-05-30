#include "PoissonSolverBase.h"
#include "PoissonSolverAMR.h"
#ifdef GPU_POISSON
#include "PoissonSolverExp.h"
#endif
#include "../SimulationData.h"

namespace cubismup3d {

std::shared_ptr<PoissonSolverBase> makePoissonSolver(SimulationData& s)
{
  if (s.poissonSolver == "iterative") 
  {
    return std::make_shared<PoissonSolverAMR>(s);
  } 
  else if (s.poissonSolver == "cuda_iterative") 
  {
#ifdef GPU_POISSON
    if (! _DOUBLE_PRECISION_ )
      throw std::runtime_error( 
          "Poisson solver: \"" + s.poissonSolver + "\" must be compiled with in double precision mode!" );
    return std::make_shared<PoissonSolverExp>(s);
#else
    throw std::runtime_error(
        "Poisson solver: \"" + s.poissonSolver + "\" must be compiled with the -DGPU_POISSON flag!"); 
#endif
  } 
  else {
    throw std::invalid_argument(
        "Poisson solver: \"" + s.poissonSolver + "\" unrecognized!");
  }
}
} // namespace cubismup3d
