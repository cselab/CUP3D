#include "Common.h"
#include "../Simulation.h"
#include "../obstacles/Obstacle.h"
#include "../obstacles/Sphere.h"

#include <mpi.h>

using namespace cubismup3d;
using namespace cubismup3d::pybindings;
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

static void bindSimulationData(py::module &m)
{
  py::class_<SimulationData>(m, "SimulationData")
      .def_readonly("CFL", &SimulationData::CFL)
      .def_readonly("BCx_flag", &SimulationData::BCx_flag, "Boundary condition in x-axis.")
      .def_readonly("BCy_flag", &SimulationData::BCy_flag, "Boundary condition in y-axis.")
      .def_readonly("BCz_flag", &SimulationData::BCz_flag, "Boundary condition in z-axis.")
      .def_readonly("extent", &SimulationData::extent)
      .def_readonly("uinf", &SimulationData::uinf)
      .def_readonly("nsteps", &SimulationData::nsteps);
}

PYBIND11_MODULE(libcubismup3d, m)
{
  using namespace py::literals;
  m.doc() = "CubismUP3D solver for incompressible Navier-Stokes";

  bindSimulationData(m);

  /* Simulation */
  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation")
      .def(py::init([](const std::vector<std::string> &argv) {
        return createSimulation(MPI_COMM_WORLD, argv);
      }), "argv"_a)
      .def_readonly("sim", &Simulation::sim, py::return_value_policy::reference)
      .def("run", &Simulation::run)
      .def("add_obstacle", &Simulation_addObstacle);

  bindObstacles(m);
}
