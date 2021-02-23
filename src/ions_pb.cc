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
    py::class_<Ions>(m, "Ions");

    py::class_<Ions_MonteCarlo,Ions>(m, "Ions_MonteCarlo")
        .def(py::init<const Twiss &, int>())
        .def(py::init<const Twiss &, std::string, int, int, bool, int>())
        .def("adjust_disp", &Ions_MonteCarlo::adjust_disp)
        .def("adjust_disp_inv", &Ions_MonteCarlo::adjust_disp_inv)
        .def("save_ions_sdds", &Ions_MonteCarlo::save_ions_sdds)
        .def("get_twiss", &Ions_MonteCarlo::get_twiss)
        .def("n_sample", &Ions_MonteCarlo::n_sample)
// TODO What the emittance getters ought to look like will become clear
//      once we have a consistent way of storing the current emittance
//        .def("emit", py::overload_cast<double&, double&, double&>(&Ions_MonteCarlo::emit))
//        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
//                vector<double>&, double&, double&, double&>(&Ions_MonteCarlo::emit))
        .def("create_samples", &Ions_MonteCarlo::create_samples);

    py::class_<Ions_SingleParticle,Ions>(m, "Ions_SingleParticle")
        .def(py::init<const Twiss &, int, int>())
        .def("single_particle_grid", &Ions_SingleParticle::single_particle_grid)
        .def("adjust_disp", &Ions_SingleParticle::adjust_disp)
        .def("adjust_disp_inv", &Ions_SingleParticle::adjust_disp_inv)
        .def("save_ions_sdds", &Ions_SingleParticle::save_ions_sdds)
        .def("get_twiss", &Ions_SingleParticle::get_twiss)
        .def("n_sample", &Ions_SingleParticle::n_sample)
//        .def("emit", py::overload_cast<double&, double&, double&>(&Ions_SingleParticle::emit))
//        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
//                vector<double>&, double&, double&, double&>(&Ions_SingleParticle::emit))
        .def("create_samples", &Ions_SingleParticle::create_samples);
}
