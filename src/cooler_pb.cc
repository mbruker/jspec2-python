#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/cooler.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_cooler(py::module& m) {
    py::class_<Cooler>(m, "Cooler")
        .def(py::init<double, double, double, double, double, double, double, double, double, double, double>(),
	    py::arg("length"),
	    py::arg("section_number"),
	    py::arg("magnetic_field"),
	    py::arg("beta_h"),
	    py::arg("beta_v"),
	    py::arg("disp_h") = 0,
	    py::arg("disp_v") = 0,
	    py::arg("alpha_h") = 0,
	    py::arg("alpha_v") = 0,
	    py::arg("der_disp_h") = 0,
	    py::arg("der_disp_v") = 0
	)
        .def("length", &Cooler::length)
        .def("section_number", &Cooler::section_number)
        .def("magnetic_filed", &Cooler::magnetic_field)
        .def("beta_h", &Cooler::beta_h)
        .def("beta_v", &Cooler::beta_v)
        .def("alpha_h", &Cooler::alpha_h)
        .def("alpha_v", &Cooler::alpha_v)
        .def("disp_h", &Cooler::disp_h)
        .def("disp_v", &Cooler::disp_v)
        .def("der_disp_h", &Cooler::der_disp_h)
        .def("der_disp_v", &Cooler::der_disp_v);
}
