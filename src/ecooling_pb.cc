#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/cooler.h"
#include "jspec2/ring.h"
#include "jspec2/ions.h"
#include "jspec2/force.h"
#include "jspec2/ecooling.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_ecooling(py::module& m) {
//    py::enum_<IonSample>(m, "IonSample", py::arithmetic())
//        .value("SINGLE_PARTICLE", IonSample::SINGLE_PARTICLE)
//        .value("MONTE_CARLO", IonSample::MONTE_CARLO);
//
//    py::enum_<ECoolRateScratch>(m, "ECoolRateScratch", py::arithmetic())
//        .value("XP_BET", ECoolRateScratch::XP_BET)
//        .value("YP_BET", ECoolRateScratch::YP_BET)
//        .value("XP", ECoolRateScratch::XP)
//        .value("YP", ECoolRateScratch::YP)
//        .value("DP_P", ECoolRateScratch::DP_P)
//        .value("V_TR", ECoolRateScratch::V_TR)
//        .value("V_LONG", ECoolRateScratch::V_LONG)
//        .value("FORCE_X", ECoolRateScratch::FORCE_X)
//        .value("FORCE_Y", ECoolRateScratch::FORCE_Y)
//        .value("FORCE_Z", ECoolRateScratch::FORCE_Z);

    py::class_<ECoolRate>(m, "ECoolRate")
        .def(py::init<>())
        .def("set_dual_force_solver", &ECoolRate::set_dual_force_solver)
        .def("set_second_force_solver", &ECoolRate::set_second_force_solver)
        .def("adjust_rate", &ECoolRate::adjust_rate)
        .def("scratch", &ECoolRate::scratch)
        .def("t_cooler", &ECoolRate::t_cooler)
        .def("set_n_long_sample", &ECoolRate::set_n_long_sample)
        .def("ecool_rate", &ECoolRate::ecool_rate);
}
