#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "jspec2/cooler.h"
#include "jspec2/ring.h"
#include "jspec2/ion_beam.h"
#include "jspec2/force.h"
#include "jspec2/ecooling.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_ecooling(py::module& m) {
    py::class_<ECoolRate>(m, "ECoolRate")
        .def(py::init<>())
        .def("set_dual_force_solver", &ECoolRate::set_dual_force_solver)
        .def("set_second_force_solver", &ECoolRate::set_second_force_solver)
        .def("adjust_rate", &ECoolRate::adjust_rate)
        .def("t_cooler", &ECoolRate::t_cooler)
        .def("set_n_long_sample", &ECoolRate::set_n_long_sample)
        .def("ecool_rate", &ECoolRate::ecool_rate);
}
