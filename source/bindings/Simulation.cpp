#include "Common.h"
#include "../Simulation.h"
#include "../obstacles/Obstacle.h"
#include "../obstacles/Sphere.h"
#include "../operators/Operator.h"

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

/// Checkpoint listener that stops the simulation if Ctrl-C was pressed.
struct SIGINTHandler : Operator
{
  using Operator::Operator;
  void operator()(double /* dt */) override
  {
    // https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-properly-handle-ctrl-c-in-long-running-functions
    if (PyErr_CheckSignals() != 0)
      throw py::error_already_set();
  }

  std::string getName() override
  {
    return "SIGINTHandler";
  }
};


static std::shared_ptr<Simulation> pyCreateSimulation(
    const std::vector<std::string> &argv,
    uintptr_t commPtr)
{
  // https://stackoverflow.com/questions/49259704/pybind11-possible-to-use-mpi4py
  // In Python, pass `MPI._addressof(comm)` as the value of the `comm` argument.
  MPI_Comm comm = commPtr ? *(MPI_Comm *)commPtr : MPI_COMM_WORLD;
  auto sim = createSimulation(comm, argv);
  sim->sim.pipeline.push_back(std::make_shared<SIGINTHandler>(sim->sim));
  return sim;
}

PYBIND11_MODULE(libcubismup3d, m)
{
  using namespace py::literals;
  m.doc() = "CubismUP3D solver for incompressible Navier-Stokes";

  bindSimulationData(m);

  /* Simulation */
  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation")
      .def(py::init(&pyCreateSimulation), "argv"_a, "comm"_a = 0)
      .def_readonly("sim", &Simulation::sim, py::return_value_policy::reference_internal)
      .def("run", &Simulation::run)
      .def("add_obstacle", &Simulation_addObstacle)
      .def("add_obstacle", &Simulation_parseAndAddObstacle);

  bindObstacles(m);
}
