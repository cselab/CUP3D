#include "Common.h"
#include "../operators/Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

namespace py = pybind11;
using namespace py::literals;

class PyOperator : public Operator
{
public:
  PyOperator(SimulationData& s, std::string name) :
    Operator{s}, name_{std::move(name)}
  { }

  void operator()(const Real dt) override
  {
    // If this fails, store the object in a permanent variable somewhere. See
    // https://github.com/pybind/pybind11/issues/1546
    // https://github.com/pybind/pybind11/pull/2839
    PYBIND11_OVERRIDE_PURE_NAME(void, Operator, "__call__", operator(), dt);
  }

  std::string getName() override
  {
    return name_;
  }

private:
  std::string name_;
};

void bindOperators(py::module &m)
{
  using namespace py::literals;
  py::class_<Operator, PyOperator, std::shared_ptr<Operator>>(m, "Operator")
    .def(py::init<SimulationData&, std::string>(), "data"_a, "name"_a)
    .def("__str__", &Operator::getName)
    .def("__repr__", &Operator::getName)
    .def_property_readonly("data", [](Operator *op) { return &op->sim; },
                           py::return_value_policy::reference_internal)
    .def("__call__", &Operator::operator(), "dt"_a);
}

CubismUP_3D_NAMESPACE_END
