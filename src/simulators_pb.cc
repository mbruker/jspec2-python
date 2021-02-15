#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <functional>
#include <tuple>
#include <vector>
#include <string>
#include <iostream>

#include "jspec2/simulator.h"
#include "jspec2/cooler.h"
#include "jspec2/ions.h"
#include "jspec2/particle_model.h"
#include "jspec2/rms_dynamic.h"
#include "jspec2/turn_by_turn.h"
#include "jspec2/force.h"
#include "jspec2/ibs.h"
#include "jspec2/ecooling.h"
#include "jspec2/luminosity.h"
#include "jspec2/ring.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

class PySink : public DataSink {
public:
    
    py::object queue;
    py::function put;
    PySink(py::object _queue);

    virtual void output_simulator_state(const Simulator::State &state) override;
};

PySink::PySink(py::object _queue)
    : queue(_queue)
{
    put = queue.attr("send");
}

void PySink::output_simulator_state(const Simulator::State &state)
{
    py::array_t<Simulator::State> state_array(1);
    *(static_cast<Simulator::State *>(state_array.request().ptr)) = state;
    std::cout << "putting t=" << state.t << std::endl;
    put(state_array);
}

void init_simulators(py::module& m) {
    py::class_<DataSink>(m, "DataSink");
    py::class_<PySink, DataSink>(m, "PySink")
        .def(py::init<py::object>());
        
    py::class_<Simulator>(m, "Simulator")
        .def("set_datasink", &Simulator::set_datasink);

    PYBIND11_NUMPY_DTYPE(Simulator::State, t, ex, ey, dp_p, sigma_s, rx, ry, rs, rx_ibs, ry_ibs, rs_ibs, rx_ecool, ry_ecool, rs_ecool, rf_voltage, luminosity);
    
        
    py::class_<ParticleModel, Simulator>(m, "ParticleModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, FrictionForceSolver *, LuminositySolver *>())
//        .def("set_ion_save", &ParticleModel::set_ion_save)
//        .def("set_output_file", &ParticleModel::set_output_file)
//        .def("set_output_intvl", &ParticleModel::set_output_intvl)
        .def("set_fixed_bunch_length", &ParticleModel::set_fixed_bunch_length)
//        .def("set_ini_time", &ParticleModel::set_ini_time)
//        .def("set_reset_time", &ParticleModel::set_reset_time)
//        .def("set_overwrite", &ParticleModel::set_overwrite)
//        .def("set_calc_lum", &ParticleModel::set_calc_lum)
        .def("run", &ParticleModel::run);

    py::class_<RMSModel, Simulator>(m, "RMSModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, FrictionForceSolver *, LuminositySolver *>())
//        .def("set_ion_save", &RMSModel::set_ion_save)
//        .def("set_output_file", &RMSModel::set_output_file)
//        .def("set_output_intvl", &RMSModel::set_output_intvl)
        .def("set_fixed_bunch_length", &RMSModel::set_fixed_bunch_length)
//        .def("set_ini_time", &RMSModel::set_ini_time)
//        .def("set_reset_time", &RMSModel::set_reset_time)
//        .def("set_overwrite", &RMSModel::set_overwrite)
//        .def("set_calc_lum", &RMSModel::set_calc_lum)
        .def("run", &RMSModel::run);

    py::class_<TurnByTurnModel, Simulator>(m, "TurnByTurnModel")
        .def(py::init<Ring &, Cooler &, IBSSolver *, ECoolRate *, FrictionForceSolver *, LuminositySolver *>())
//        .def("set_ion_save", &TurnByTurnModel::set_ion_save)
//        .def("set_output_file", &TurnByTurnModel::set_output_file)
//        .def("set_output_intvl", &TurnByTurnModel::set_output_intvl)
        .def("set_fixed_bunch_length", &TurnByTurnModel::set_fixed_bunch_length)
//        .def("set_ini_time", &TurnByTurnModel::set_ini_time)
//        .def("set_reset_time", &TurnByTurnModel::set_reset_time)
//        .def("set_overwrite", &TurnByTurnModel::set_overwrite)
//        .def("set_calc_lum", &TurnByTurnModel::set_calc_lum)
        .def("run", &TurnByTurnModel::run);
}
