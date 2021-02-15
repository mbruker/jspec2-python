#include <pybind11/pybind11.h>
#include "jspec2/cooler.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_cooler(py::module& m) {
    py::class_<Cooler>(m, "Cooler")
        .def(py::init<double, double, double, const Twiss &, double>(),
	    py::arg("length"),
	    py::arg("section_number") = 1,
	    py::arg("magnetic_field"),
	    py::arg("twiss"),
        py::arg("pipe_radius") = 0
	)
        .def_property_readonly("length", &Cooler::length)
        .def_property_readonly("section_number", &Cooler::section_number)
        .def_property_readonly("magnetic_field", &Cooler::magnetic_field)
        .def_property_readonly("twiss", &Cooler::twiss)
        .def_property_readonly("pipe_radius", &Cooler::pipe_radius);
}
