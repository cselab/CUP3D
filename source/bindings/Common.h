#pragma once

#include "../Base.h"
#include <array>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

CubismUP_3D_NAMESPACE_BEGIN

class Obstacle;
class Simulation;

void bindSimulationData(pybind11::module &m);
void bindSimulation(pybind11::module &m);
void bindObstacles(pybind11::module &m);
void bindOperators(pybind11::module &m);
void pySimulationAddObstacle(Simulation &s, std::shared_ptr<Obstacle>);
void pySimulationParseAndAddObstacle(Simulation &S, pybind11::object obstacle_args);

CubismUP_3D_NAMESPACE_END
