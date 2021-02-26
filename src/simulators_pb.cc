#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "jspec2/simulator.h"
#include "jspec2/cooler.h"
#include "jspec2/ion_beam.h"
#include "jspec2/particle_model.h"
#include "jspec2/rms_dynamic.h"
#include "jspec2/turn_by_turn.h"
#include "jspec2/ibs.h"
#include "jspec2/ecooling.h"
#include "jspec2/luminosity.h"
#include "jspec2/ring.h"
#include "jspec2/datasink.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_simulators(py::module& m)
{
    py::class_<Simulator>(m, "Simulator")
        .def("set_datasink", &Simulator::set_datasink);

    py::class_<ParticleModel, Simulator>(m, "ParticleModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, LuminositySolver *>())
        .def("set_fixed_bunch_length", &ParticleModel::set_fixed_bunch_length)
//        .def("set_ini_time", &ParticleModel::set_ini_time)
//        .def("set_reset_time", &ParticleModel::set_reset_time)
        .def("run", &ParticleModel::run);

    py::class_<RMSModel, Simulator>(m, "RMSModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, LuminositySolver *>())
        .def("set_fixed_bunch_length", &RMSModel::set_fixed_bunch_length)
//        .def("set_ini_time", &RMSModel::set_ini_time)
//        .def("set_reset_time", &RMSModel::set_reset_time)
        .def("run", &RMSModel::run);

    py::class_<TurnByTurnModel, Simulator>(m, "TurnByTurnModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, LuminositySolver *>())
        .def("set_fixed_bunch_length", &TurnByTurnModel::set_fixed_bunch_length)
//        .def("set_ini_time", &TurnByTurnModel::set_ini_time)
//        .def("set_reset_time", &TurnByTurnModel::set_reset_time)
        .def("run", &TurnByTurnModel::run);
}
