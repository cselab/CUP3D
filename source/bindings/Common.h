#ifndef CUBIMSUP3D_BINDINGS_COMMON_H
#define CUBIMSUP3D_BINDINGS_COMMON_H

#include <array>
#include <pybind11/stl.h>

namespace cubismup3d {
namespace pytypes {

using int3 = std::array<int, 3>;
using double3 = std::array<double, 3>;

}  // namespace pytypes
}  // namespace cubismup3d

#endif
