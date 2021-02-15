#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>
#include <string>


#include "jspec2/cooler.h"
#include "jspec2/ions.h"
#include "jspec2/ring.h"
#include "jspec2/twiss.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_ions(py::module& m) {
    py::enum_<IonSampleType>(m, "IonSampleType", py::arithmetic())
        .value("MONTE_CARLO", IonSampleType::MONTE_CARLO)
        .value("USER_DEFINE", IonSampleType::USER_DEFINE)
        .value("SINGLE_PARTICLE", IonSampleType::SINGLE_PARTICLE);

    py::enum_<Phase>(m, "Phase", py::arithmetic())
        .value("X_BET", Phase::X_BET)
        .value("XP_BET", Phase::XP_BET)
        .value("Y_BET", Phase::Y_BET)
        .value("YP_BET", Phase::YP_BET)
        .value("X", Phase::X)
        .value("Y", Phase::Y)
        .value("XP", Phase::XP)
        .value("YP", Phase::YP)
        .value("DS", Phase::DS)
        .value("DP_P", Phase::DP_P);

    py::class_<Ions>(m, "Ions");
//        .def(py::init<const Twiss &>());

    py::class_<Ions_MonteCarlo,Ions>(m, "Ions_MonteCarlo")
        .def(py::init<const Twiss &, int>())
        .def(py::init<const Twiss &, std::string, int, int, bool, int>())
        .def("adjust_disp", &Ions_MonteCarlo::adjust_disp)
        .def("adjust_disp_inv", &Ions_MonteCarlo::adjust_disp_inv)
        .def("save_ions_sdds", &Ions_MonteCarlo::save_ions_sdds)
        .def("get_cdnt", &Ions_MonteCarlo::get_cdnt)
        .def("cdnt", &Ions_MonteCarlo::cdnt)
        .def("get_twiss", &Ions_MonteCarlo::get_twiss)
        .def("n_sample", &Ions_MonteCarlo::n_sample)
        .def("sample_type", &Ions_MonteCarlo::sample_type)
        .def("center", &Ions_MonteCarlo::center)
        .def("emit", py::overload_cast<double&, double&, double&>(&Ions_MonteCarlo::emit))
        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
                vector<double>&, double&, double&, double&>(&Ions_MonteCarlo::emit))
        .def("create_samples", &Ions_MonteCarlo::create_samples);

    py::class_<Ions_SingleParticle,Ions>(m, "Ions_SingleParticle")
        .def(py::init<const Twiss &, int, int>())
        .def("single_particle_grid", &Ions_SingleParticle::single_particle_grid)
        .def("adjust_disp", &Ions_SingleParticle::adjust_disp)
        .def("adjust_disp_inv", &Ions_SingleParticle::adjust_disp_inv)
        .def("save_ions_sdds", &Ions_SingleParticle::save_ions_sdds)
        .def("get_cdnt", &Ions_SingleParticle::get_cdnt)
        .def("cdnt", &Ions_SingleParticle::cdnt)
        .def("get_twiss", &Ions_SingleParticle::get_twiss)
        .def("n_sample", &Ions_SingleParticle::n_sample)
        .def("sample_type", &Ions_SingleParticle::sample_type)
        .def("center", &Ions_SingleParticle::center)
        .def("emit", py::overload_cast<double&, double&, double&>(&Ions_SingleParticle::emit))
        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
                vector<double>&, double&, double&, double&>(&Ions_SingleParticle::emit))
        .def("create_samples", &Ions_SingleParticle::create_samples);

    m.def("emit_p", &emit_p);
    m.def("emit", &emit);
    m.def("adjust_disp", &adjust_disp);
    m.def("adjust_disp_inv", &adjust_disp_inv);
}
