#include "Common.h"
#include "../Simulation.h"
#include "../operators/ImportExportUniform.h"

#include <pybind11/numpy.h>

CubismUP_3D_NAMESPACE_BEGIN

namespace py = pybind11;
using namespace py::literals;

namespace {

/// View into one of the 8 fluid components.
struct PyFieldView
{
  Simulation *sim;
  int component;

  py::str str() const
  {
    return py::str("_FieldView(sim={}, component={})").format(sim, component);
  }

  void loadUniform(py::array_t<Real, py::array::c_style> in)
  {
    _checkShape(in);
    importGridFieldFromUniformMatrix(sim->sim.grid, component, in.data());
  }

  auto toUniform(py::array_t<Real, py::array::c_style> out) const
  {
    _checkShape(out);
    exportGridFieldToUniformMatrix(sim->sim.grid, component, out.mutable_data());
    return out;
  }

  void _checkShape(const py::array_t<Real, py::array::c_style> &out) const
  {
    const auto numCells = sim->sim.grid->getMaxMostRefinedCells();
    if (out.ndim() != 3
        || out.shape(2) != numCells[0]
        || out.shape(1) != numCells[1]
        || out.shape(0) != numCells[2]) {
      throw py::type_error(
          py::str("expected shape ({}, {}, {}), got {}").format(
            numCells[2], numCells[1], numCells[0],
            std::vector(out.shape(), out.shape() + out.ndim())));
    }
  }
};

}  // anonymous namespace

static auto fieldViewProperty(int which)
{
  return [which](PyFieldsView view)
  {
    return PyFieldView{view.sim, which};
  };
}

void bindFields(py::module &m)
{
  py::class_<PyFieldView>(m, "_FieldView")
    .def("__str__", &PyFieldView::str)
    .def("__repr__", &PyFieldView::str)
    .def("load_uniform", &PyFieldView::loadUniform, "field"_a.noconvert(),
         "Import field from a contiguous matrix representing the fully "
         "refined uniform grid. This function does not affect the AMR mesh "
         "itself.")
    .def("to_uniform", &PyFieldView::toUniform, "out"_a.noconvert(),
         "Copy and interpolate the locally available part of the AMR field "
         "into a uniform matrix. Cells belonging to other ranks are left "
         "unchanged. The `out` array must be contiguous in the C-order.");

  py::class_<PyFieldsView>(m, "_FieldsView")
    .def_property_readonly("chi", fieldViewProperty(FE_CHI))
    .def_property_readonly("u", fieldViewProperty(FE_U))
    .def_property_readonly("v", fieldViewProperty(FE_V))
    .def_property_readonly("w", fieldViewProperty(FE_W))
    .def_property_readonly("p", fieldViewProperty(FE_P))
    .def_property_readonly("tmpU", fieldViewProperty(FE_TMPU))
    .def_property_readonly("tmpV", fieldViewProperty(FE_TMPV))
    .def_property_readonly("tmpW", fieldViewProperty(FE_TMPW));
}

CubismUP_3D_NAMESPACE_END

