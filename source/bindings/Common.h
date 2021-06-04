#ifndef CUBIMSUP3D_BINDINGS_COMMON_H
#define CUBIMSUP3D_BINDINGS_COMMON_H

#include "../Base.h"
#include "../obstacles/Obstacle.h"

#include <array>
#include <pybind11/stl.h>

CubismUP_3D_NAMESPACE_BEGIN

class Simulation;

namespace pybindings {

void bindObstacles(pybind11::module &m);
void Simulation_addObstacle(Simulation &s, std::shared_ptr<Obstacle>);
void Simulation_parseAndAddObstacle(Simulation &S, pybind11::object obstacle_args);

}  // namespace pybindings

CubismUP_3D_NAMESPACE_END

#endif
