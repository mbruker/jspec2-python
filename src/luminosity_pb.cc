#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>
#include <string>

#include "jspec2/luminosity.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_luminosity(py::module& m) {
     py::class_<CollidingBeam>(m, "CollidingBeam")
        .def(py::init<>())
        .def_readwrite("np",&CollidingBeam::np)
        .def_readwrite("sigma_x",&CollidingBeam::sigma_x)
        .def_readwrite("sigma_y",&CollidingBeam::sigma_y)
        .def_readwrite("geo_emit_x",&CollidingBeam::geo_emit_x)
        .def_readwrite("geo_emit_y",&CollidingBeam::geo_emit_y)
        .def_readwrite("bet_x_star",&CollidingBeam::bet_x_star)
        .def_readwrite("bet_y_star",&CollidingBeam::bet_y_star)
        .def_readwrite("bet_x_max",&CollidingBeam::bet_x_max)
        .def_readwrite("bet_y_max",&CollidingBeam::bet_y_max)
        .def_readwrite("focus_length_x",&CollidingBeam::focus_length_x)
        .def_readwrite("focus_length_y",&CollidingBeam::focus_length_y)
        .def_readwrite("aper_x",&CollidingBeam::aper_x)
        .def_readwrite("aper_y",&CollidingBeam::aper_y)
        .def_readwrite("adjust_bet",&CollidingBeam::adjust_bet);

    py::class_<LuminositySolver>(m, "LuminositySolver")
        .def("set_distance", &LuminositySolver::set_distance)
        .def("set_freq", &LuminositySolver::set_freq)
        .def("set_match", &LuminositySolver::set_match)
        .def("set_geo_emit", &LuminositySolver::set_geo_emit)
        .def("set_beam_size", &LuminositySolver::set_beam_size)
        .def("set_particle_number", &LuminositySolver::set_particle_number)
        .def("set_bet_star", &LuminositySolver::set_bet_star)
        .def("set_bet_max", &LuminositySolver::set_bet_max)
        .def("match", &LuminositySolver::match)
        .def("set_focus_length", &LuminositySolver::set_focus_length)
        .def("set_aperture", &LuminositySolver::set_aperture)
        .def("set_adjust_bet", &LuminositySolver::set_adjust_bet)
        .def("set_aper_ratio", &LuminositySolver::set_aper_ratio)
        .def("set_use_ion_emit", &LuminositySolver::set_use_ion_emit)
        .def("use_ion_emit", &LuminositySolver::use_ion_emit)
        .def("bet_x_max", &LuminositySolver::bet_x_max)
        .def("bet_y_max", &LuminositySolver::bet_y_max)
        .def("bet_x_star", &LuminositySolver::bet_x_star)
        .def("bet_y_star", &LuminositySolver::bet_y_star)
        .def("luminosity", &LuminositySolver::luminosity);

}
