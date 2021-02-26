#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "jspec2/cooler.h"
#include "jspec2/datasink.h"
#include "jspec2/ring.h"
#include "jspec2/ion_beam.h"
#include "jspec2/force.h"
#include "jspec2/ecooling.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_ecooling(py::module& m) {
    py::class_<ECoolRate>(m, "ECoolRate")
        .def(py::init<FrictionForceSolver *, FrictionForceSolver *>(),
             py::arg("force_solver"),
             py::arg("longitudinal_force_solver") = nullptr
        )
        .def("set_n_long_sample", &ECoolRate::set_n_long_sample)
        .def("set_force_datasink", &ECoolRate::set_force_datasink)
        .def("ecool_rate", &ECoolRate::ecool_rate);
        
    py::class_<ForceCurve, ECoolRate>(m, "ForceCurve")
        .def("output_force", &ForceCurve::output_force);
}
