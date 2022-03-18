#include "Common.h"

#include <mpi.h>

using namespace cubismup3d;
using namespace pybind11::literals;
namespace py = pybind11;

namespace {

/* Ensure that we load highest thread level we need. */
struct CUPMPILoader
{
  CUPMPILoader()
  {
    int flag, provided;
    MPI_Initialized(&flag);
    if (!flag) {
      MPI_Init_thread(0, nullptr, MPI_THREAD_MULTIPLE, &provided);
    } else {
      MPI_Query_thread(&provided);
    }
#ifdef CUP_ASYNC_DUMP
    const auto SECURITY = MPI_THREAD_MULTIPLE;
#else
    const auto SECURITY = MPI_THREAD_FUNNELED;
#endif
    if (provided >= SECURITY)
      return;
    if (!flag)
      fprintf(stderr, "Error: MPI implementation does not have the required thread support!\n");
    else
      fprintf(stderr, "Error: MPI does not implement or not initialized with the required thread support!\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
} cup_mpi_loader;

}  // anonymous namespace

PYBIND11_MODULE(libcubismup3d, m)
{
  m.doc() = "CubismUP3D solver for incompressible Navier-Stokes";

  bindFields(m);
  bindSimulationData(m);
  bindSimulation(m);
  bindObstacles(m);
  bindOperators(m);
}
