#include <mpi.h>
#include "bindings/Common.h"
#include "Simulation.h"

using namespace cubismup3d::pytypes;
using namespace pybind11::literals;

namespace {

struct CUP_MPI_Loader {
  CUP_MPI_Loader() {
    int flag, provided;
    MPI_Initialized(&flag);
    if (!flag) {
      MPI_Init_thread(0, nullptr, MPI_THREAD_MULTIPLE, &provided);
    } else {
      MPI_Query_thread(&provided);
    }
#ifdef DUMPGRID
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

}  // namespace

namespace py = pybind11;
PYBIND11_MODULE(_cubismup3d, m) {

    m.doc() = "CubismUP3D solver for incompressible Navier-Stokes";

    // This should support everything that ArgumentParser version of the
    // Simulation constructor has.
    py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation")
        .def(py::init<int3, int3,
                      MPI_Comm,
                      int, double,
                      double, double, double, double,
                      double3,
                      bool, bool,
                      bool, bool,
#ifndef _UNBOUNDED_FFT_
                      double /* fadeOutLength */,
#endif
                      int, double,
                      const std::string &,
                      bool>(),
             py::return_value_policy::take_ownership,
             "cells"_a, "nproc"_a = int3{-1, -1, -1},
             "comm"_a = MPI_COMM_WORLD,
             "nsteps"_a = 0, "endTime"_a = 0.0,
             "nu"_a = 0.0, "CFL"_a = 0.1, "lambda_"_a = 0.0, "DLM"_a = 1.0,
             "uinf"_a = double3{0.0, 0.0, 0.0},
             "verbose"_a = true, "computeDissipation"_a = false,
             "dump3D"_a = true, "dump2D"_a = false,
#ifndef _UNBOUNDED_FFT_
             "fadeOutLength"_a = 0.005,
#endif
             "saveFreq"_a = 0, "saveTime"_a = 0.0,
             "path4serialization"_a = std::string("./"),
             "restart"_a = false, R"(
            Documentation....
        )")
        .def_readonly("bpdx", &Simulation::bpdx, "Blocks in x-direction.")
        .def_readonly("bpdy", &Simulation::bpdy, "Blocks in y-direction.")
        .def_readonly("bpdz", &Simulation::bpdz, "Blocks in z-direction.")
        .def_readonly("nprocsx", &Simulation::nprocsx)
        .def_readonly("nprocsy", &Simulation::nprocsy)
        .def_readonly("nprocsz", &Simulation::nprocsz);
}
