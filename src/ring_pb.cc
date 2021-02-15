#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/ring.h"
#include "jspec2/beam.h"

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
        .def(py::init<const Lattice&, const Beam *>())
        .def_readwrite("tunes", &Ring::tunes)
        .def_readwrite("rf", &Ring::rf)
        .def("beta_s", &Ring::beta_s)
        .def("circ", &Ring::circ)
        .def("f0", &Ring::f0)
        .def("w0", &Ring::w0)
        .def("slip_factor", &Ring::slip_factor)
        .def("calc_rf_voltage", &Ring::calc_rf_voltage)
        .def("calc_sync_tune_by_rf", &Ring::calc_sync_tune_by_rf)
        .def("beam", &Ring::beam)
        .def("lattice", &Ring::lattice)
        .def("set_rf", &Ring::set_rf);

    py::class_<Tunes>(m, "Tunes")
        .def(py::init<>())
        .def_readwrite("qx", &Tunes::qx)
        .def_readwrite("qy", &Tunes::qy)
        .def_readwrite("qs", &Tunes::qs);

    py::class_<RF>(m, "RF")
        .def(py::init<>())
        .def_readwrite("v", &RF::v)
        .def_readwrite("h", &RF::h)
        .def_readwrite("phi", &RF::phi)
        .def_readwrite("gamma_tr", &RF::gamma_tr);

    py::class_<Twiss>(m, "Twiss")
        .def(py::init<>())
        .def_readwrite("bet_x", &Twiss::bet_x)
        .def_readwrite("alf_x", &Twiss::alf_x)
        .def_readwrite("disp_x", &Twiss::disp_x)
        .def_readwrite("disp_dx", &Twiss::disp_dx)
        .def_readwrite("bet_y", &Twiss::bet_y)
        .def_readwrite("alf_y", &Twiss::alf_y)
        .def_readwrite("disp_y", &Twiss::disp_y)
        .def_readwrite("disp_dy", &Twiss::disp_dy);

}

