#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>
#include <string>


#include "jspec2/cooler.h"
#include "jspec2/ion_beam.h"
#include "jspec2/ring.h"
#include "jspec2/twiss.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_ions(py::module& m) {
    py::class_<IonBeam>(m, "IonBeam")
        .def_property_readonly("charge_number", &IonBeam::charge_number)
        .def_property_readonly("mass", &IonBeam::mass)
        .def_property_readonly("kinetic_energy", &IonBeam::kinetic_energy)
        .def_property_readonly("beta", &IonBeam::beta)
        .def_property_readonly("gamma", &IonBeam::gamma);

    py::class_<IonBeam_MonteCarlo,IonBeam>(m, "IonBeam_MonteCarlo")
        .def(py::init<const Twiss &, int, double, double, double, double, double, double, double, int>(),
             py::arg("twiss"),
             py::arg("charge_number"),
             py::arg("mass"),
             py::arg("kinetic_energy"),
             py::arg("emit_nx"),
             py::arg("emit_ny"),
             py::arg("dp_p"),
             py::arg("sigma_s") = 0,
             py::arg("n_particle"),
             py::arg("n_sample")
            )
        .def(py::init<const Twiss &, int, double, double, double, double, double, double, double, std::string, int, int, bool, int>(),
             py::arg("twiss"),
             py::arg("charge_number"),
             py::arg("mass"),
             py::arg("kinetic_energy"),
             py::arg("emit_nx"),
             py::arg("emit_ny"),
             py::arg("dp_p"),
             py::arg("sigma_s") = 0,
             py::arg("n_particle"),
             py::arg("filename"),
             py::arg("n_sample"),
             py::arg("skip") = 0,
             py::arg("binary") = false,
             py::arg("n_buffer") = 1000
            )
        .def("adjust_disp", &IonBeam_MonteCarlo::adjust_disp)
        .def("adjust_disp_inv", &IonBeam_MonteCarlo::adjust_disp_inv)
        .def("save_ions_sdds", &IonBeam_MonteCarlo::save_ions_sdds)
        .def("get_twiss", &IonBeam_MonteCarlo::get_twiss)
        .def("n_sample", &IonBeam_MonteCarlo::n_sample)
// TODO What the emittance getters ought to look like will become clear
//      once we have a consistent way of storing the current emittance
//        .def("emit", py::overload_cast<double&, double&, double&>(&IonBeam_MonteCarlo::emit))
//        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
//                vector<double>&, double&, double&, double&>(&IonBeam_MonteCarlo::emit))
        .def("create_samples", &IonBeam_MonteCarlo::create_samples);

    py::class_<IonBeam_SingleParticle,IonBeam>(m, "IonBeam_SingleParticle")
        .def(py::init<const Twiss &, int, double, double, double, double, double, double, double, int, int>(),
             py::arg("twiss"),
             py::arg("charge_number"),
             py::arg("mass"),
             py::arg("kinetic_energy"),
             py::arg("emit_nx"),
             py::arg("emit_ny"),
             py::arg("dp_p"),
             py::arg("sigma_s") = 0,
             py::arg("n_particle"),
             py::arg("n_tr"),
             py::arg("n_l")
            )
        .def("single_particle_grid", &IonBeam_SingleParticle::single_particle_grid)
        .def("adjust_disp", &IonBeam_SingleParticle::adjust_disp)
        .def("adjust_disp_inv", &IonBeam_SingleParticle::adjust_disp_inv)
        .def("save_ions_sdds", &IonBeam_SingleParticle::save_ions_sdds)
        .def("get_twiss", &IonBeam_SingleParticle::get_twiss)
        .def("n_sample", &IonBeam_SingleParticle::n_sample)
//        .def("emit", py::overload_cast<double&, double&, double&>(&IonBeam_SingleParticle::emit))
//        .def("emit", py::overload_cast<vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&,
//                vector<double>&, double&, double&, double&>(&IonBeam_SingleParticle::emit))
        .def("create_samples", &IonBeam_SingleParticle::create_samples);
}
