#include "Common.h"
#include "../Obstacles/Obstacle.h"
#include "../Simulation.h"
#include "../SimulationData.h"
#include "../operators/Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

namespace py = pybind11;
using namespace py::literals;

namespace {

/// Checkpoint listener that stops the simulation if Ctrl-C was pressed.
struct SIGINTHandler : Operator
{
  using Operator::Operator;
  void operator()(Real /* dt */) override
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

}  // anonymous namespace


void bindSimulationData(py::module &m)
{
  py::class_<SimulationData>(m, "SimulationData")
    .def_readonly("CFL", &SimulationData::CFL)
    .def_readonly("BCx_flag", &SimulationData::BCx_flag, "boundary condition in x-axis")
    .def_readonly("BCy_flag", &SimulationData::BCy_flag, "boundary condition in y-axis")
    .def_readonly("BCz_flag", &SimulationData::BCz_flag, "boundary condition in z-axis")
    .def_readonly("extents", &SimulationData::extents, "domain size")
    .def_readonly("uinf", &SimulationData::uinf, "frame of reference velocity")
    .def_readonly("nu", &SimulationData::nu, "kinematic viscosity")
    .def_readwrite("nsteps", &SimulationData::nsteps, "maximum number of steps")
    .def_readwrite("step", &SimulationData::step, "current step")
    .def_readwrite("time", &SimulationData::time, "current time");
}

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

void bindSimulation(py::module &m)
{
  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation")
    .def(py::init(&pyCreateSimulation), "argv"_a, "comm"_a = 0)
    .def_readonly("data", &Simulation::sim,
                  py::return_value_policy::reference_internal)
    .def_property_readonly("fields", [](Simulation *sim) { return PyFieldsView{sim}; })
    .def_property_readonly("obstacles", &Simulation::getObstacleVector)
    .def("add_obstacle", &pySimulationAddObstacle)
    .def("add_obstacle", &pySimulationParseAndAddObstacle)
    .def("adapt_mesh", &Simulation::adaptMesh)
    .def("compute_vorticity", &Simulation::computeVorticity,
         "compute vorticity and store to tmpU, tmpV and tmpW fields")
    .def("insert_operator", &Simulation::insertOperator, "op"_a)
    .def("run", &Simulation::run);
}

CubismUP_3D_NAMESPACE_END
