#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/ring.h"
#include "jspec2/twiss.h"
#include "jspec2/ion_beam.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_ring(py::module& m ){
    py::class_<Lattice>(m, "Lattice")
        .def(py::init<std::string &>())
        .def("s", (double(Lattice::*)(int)) &Lattice::s, "i"_a)
        .def("betx", (double(Lattice::*)(int)) &Lattice::betx, "i"_a)
        .def("alfx", (double(Lattice::*)(int)) &Lattice::alfx, "i"_a)
        .def("mux", (double(Lattice::*)(int)) &Lattice::mux, "i"_a)
        .def("dx", (double(Lattice::*)(int)) &Lattice::dx, "i"_a)
        .def("dpx", (double(Lattice::*)(int)) &Lattice::dpx, "i"_a)
        .def("bety", (double(Lattice::*)(int)) &Lattice::bety, "i"_a)
        .def("alfy", (double(Lattice::*)(int)) &Lattice::alfy, "i"_a)
        .def("muy", (double(Lattice::*)(int)) &Lattice::muy, "i"_a)
        .def("dy", (double(Lattice::*)(int)) &Lattice::dy, "i"_a)
        .def("dpy", (double(Lattice::*)(int)) &Lattice::dpy, "i"_a)
        .def("n_element", &Lattice::n_element)
        .def("l_element", (double(Lattice::*)(int)) &Lattice::l_element, "i"_a)
        .def("circ", &Lattice::circ);

    py::class_<Ring>(m, "Ring")
        .def(py::init<const Lattice&, const IonBeam *, double, double, double, double, int, double, double>(),
             py::arg("lattice"),
             py::arg("ion_beam"),
             py::arg("qx") = 0,
             py::arg("qy") = 0,
             py::arg("qs") = 0,
             py::arg("rf_voltage") = 0,
             py::arg("rf_h") = 1,
             py::arg("rf_phi") = 0,
             py::arg("gamma_tr") = 0
        )
        .def_property_readonly("circ", &Ring::circ)
        .def_property_readonly("f0", &Ring::f0)
        .def_property_readonly("w0", &Ring::w0)
        .def_property_readonly("slip_factor", &Ring::slip_factor)
        .def_property_readonly("qx", &Ring::qx)
        .def_property_readonly("qy", &Ring::qy)
        .def_property_readonly("qs", &Ring::qs)
        .def_property_readonly("rf_voltage", &Ring::rf_voltage)
        .def_property_readonly("rf_h", &Ring::rf_h)
        .def_property_readonly("rf_phi", &Ring::rf_phi)
        .def_property_readonly("gamma_tr", &Ring::gamma_tr);

    py::class_<Twiss>(m, "Twiss")
        .def(py::init<double, double, double, double, double, double, double, double>(),
             py::arg("beta_x"),
             py::arg("beta_y"),
             py::arg("alpha_x") = 0,
             py::arg("alpha_y") = 0,
             py::arg("disp_x") = 0,
             py::arg("disp_y") = 0,
             py::arg("disp_dx") = 0,
             py::arg("disp_dy") = 0
        )
        .def_readonly("bet_x", &Twiss::bet_x)
        .def_readonly("alf_x", &Twiss::alf_x)
        .def_readonly("disp_x", &Twiss::disp_x)
        .def_readonly("disp_dx", &Twiss::disp_dx)
        .def_readonly("bet_y", &Twiss::bet_y)
        .def_readonly("alf_y", &Twiss::alf_y)
        .def_readonly("disp_y", &Twiss::disp_y)
        .def_readonly("disp_dy", &Twiss::disp_dy);
}

