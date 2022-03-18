#include "Common.h"
#include "../Simulation.h"
#include "../obstacles/ObstacleVector.h"
#include "../obstacles/Sphere.h"
#include <pybind11/functional.h>
#include <pybind11/stl.h>

using namespace pybind11::literals;
namespace py = pybind11;

CubismUP_3D_NAMESPACE_BEGIN

#define ATTR_FROM_KWARGS(o, item) \
      do { \
        (o).item = py::cast<decltype((o).item)>(kwargs_pop(#item, (o).item)); \
      } while(0)
static ObstacleArguments popObstacleArguments(py::object &kwargs_pop)
{
  ObstacleArguments o;
  ATTR_FROM_KWARGS(o, length);
  ATTR_FROM_KWARGS(o, position);
  ATTR_FROM_KWARGS(o, quaternion);
  ATTR_FROM_KWARGS(o, enforcedVelocity);
  ATTR_FROM_KWARGS(o, bForcedInSimFrame);
  ATTR_FROM_KWARGS(o, bFixFrameOfRef);
  ATTR_FROM_KWARGS(o, bFixToPlanar);
  return o;
}
static SphereArguments popSphereArguments(py::object &kwargs_pop)
{
  SphereArguments s(py::cast<double>(kwargs_pop("radius")));
  ATTR_FROM_KWARGS(s, umax);
  ATTR_FROM_KWARGS(s, tmax);
  ATTR_FROM_KWARGS(s, accel_decel);
  ATTR_FROM_KWARGS(s, bHemi);
  return s;
}
#undef ATTR_FROM_KWARGS

static std::shared_ptr<ObstacleAndSphereArguments> createObstacleAndSphereArguments(
    double radius, double umax, double tmax, bool accel_decel, bool bHemi,
    py::kwargs kwargs)
{
  // for (const auto &item : kwargs) {
  //   fprintf(stderr, "[%s]=%s\n",
  //       py::cast<std::string>(py::str(item.first)).c_str(),
  //       py::cast<std::string>(py::str(item.second)).c_str());
  // }
  py::object kwargs_pop = kwargs.attr("pop");

  double length = py::cast<double>(kwargs_pop("length", 0.0));
  if (radius > 0 && length > 0 && radius != 0.5 * length) {
    throw std::invalid_argument("cannot specify both `radius` and `length`");
  } else if (radius <= 0 && length <= 0) {
    throw std::invalid_argument("expected a `radius` or `length`");
  } else if (radius > 0) {
    length = 2.0 * radius;
  } else {
    radius = 0.5 * length;
  }

  kwargs["radius"] = radius;
  kwargs["length"] = length;

  // TODO: throw std:invalid_argument if kwargs not empty after construction.
  return std::make_shared<ObstacleAndSphereArguments>(
      popObstacleArguments(kwargs_pop),
      popSphereArguments(kwargs_pop));
}

void bindObstacles(py::module &m)
{
  /* ObstacleArguments */
  py::class_<ObstacleArguments, std::shared_ptr<ObstacleArguments>>(m, "ObstacleArguments")
    .def_readwrite("length", &ObstacleArguments::length)
    .def_readwrite("position", &ObstacleArguments::position)
    .def_readwrite("enforcedVelocity", &ObstacleArguments::enforcedVelocity)
    .def_readwrite("bForcedInSimFrame", &ObstacleArguments::bForcedInSimFrame)
    .def_readwrite("bFixFrameOfRef", &ObstacleArguments::bFixFrameOfRef)
    .def_readwrite("bFixToPlanar", &ObstacleArguments::bFixToPlanar);

  /* Obstacle */
  py::class_<Obstacle, std::shared_ptr<Obstacle>>(m, "Obstacle");

  /* SphereArguments */
  SphereArguments sa(0.1);  // Default arguments.
  py::class_<SphereArguments, std::shared_ptr<SphereArguments>>(m, "SphereArguments")
    .def_readonly("radius", &SphereArguments::radius)
    .def_readwrite("umax", &SphereArguments::umax)
    .def_readwrite("tmax", &SphereArguments::tmax)
    .def_readwrite("accel_decel", &SphereArguments::accel_decel)
    .def_readwrite("bHemi", &SphereArguments::bHemi);

  /* ObstacleAndSphereArguments */
  py::class_<ObstacleAndSphereArguments,
             ObstacleArguments,
             SphereArguments,
             std::shared_ptr<ObstacleAndSphereArguments>>(m, "Sphere")
      .def(py::init(&createObstacleAndSphereArguments),
           "radius"_a,
           "umax"_a = sa.umax,
           "tmax"_a = sa.tmax,
           "accel_decel"_a = sa.accel_decel,
           "bHemi"_a = sa.bHemi);

  /* Sphere */
  py::class_<Sphere, Obstacle, std::shared_ptr<Sphere>>(m, "SphereObstacle")
    .def(py::init<SimulationData &, ObstacleAndSphereArguments>());
}

void pySimulationAddObstacle(Simulation &s, std::shared_ptr<Obstacle> obstacle)
{
  s.sim.obstacle_vector->addObstacle(std::move(obstacle));
}

void pySimulationParseAndAddObstacle(Simulation &S, pybind11::object obstacle_args)
{
  if (py::isinstance<ObstacleAndSphereArguments>(obstacle_args)) {
    auto args = py::cast<ObstacleAndSphereArguments>(obstacle_args);
    S.sim.obstacle_vector->addObstacle(std::make_shared<Sphere>(S.sim, args));
  } else {
    throw std::invalid_argument(py::str(obstacle_args));
  }
}

CubismUP_3D_NAMESPACE_END
