#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/electron_beam.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_electron_beam(py::module &m) {
    py::class_<ElectronBeam>(m, "ElectronBeam")
        .def("velocity", &ElectronBeam::velocity)
        .def("temperature", &ElectronBeam::temperature)
        .def("charge_number", &ElectronBeam::charge_number)
        .def("mass", &ElectronBeam::mass)
        .def("mass_number", &ElectronBeam::mass_number)
        .def("mass_SI", &ElectronBeam::mass_SI)
        .def("kinetic_energy", &ElectronBeam::kinetic_energy)
        .def("gamma", &ElectronBeam::gamma)
        .def("beta", &ElectronBeam::beta)
        .def("bunched", &ElectronBeam::bunched)
        .def("set_p_shift", &ElectronBeam::set_p_shift)
        .def("set_v_shift", &ElectronBeam::set_v_shift)
        .def("p_shift", &ElectronBeam::p_shift)
        .def("v_shift", &ElectronBeam::v_shift)
        .def("set_cv_l", &ElectronBeam::set_cv_l)
        .def("cv_l", &ElectronBeam::cv_l)
        .def("shape", &ElectronBeam::shape)
        .def("length", &ElectronBeam::length)
        .def("neutral", &ElectronBeam::neutral)
        .def("set_kinetic_energy", &ElectronBeam::set_kinetic_energy)
        .def("set_gamma", &ElectronBeam::set_gamma)
        .def("set_beta", &ElectronBeam::set_beta)
        .def("set_center", (void (ElectronBeam::*)(double, double, double)) &ElectronBeam::set_center, "cx"_a, "cy"_a, "cz"_a)
        .def("set_tpr", &ElectronBeam::set_tpr)
        .def("set_v_rms", &ElectronBeam::set_v_rms)
        .def("set_v_avg", &ElectronBeam::set_v_avg)
        .def("set_neutral", &ElectronBeam::set_neutral)
        .def("set_multi_bunches", &ElectronBeam::set_multi_bunches)
        .def("multi_bunches", &ElectronBeam::multi_bunches)
        .def_property_readonly("v_rms_l", &ElectronBeam::get_v_rms_l)
        .def_property_readonly("v_rms_tr", &ElectronBeam::get_v_rms_tr)
        .def("set_n_bunches", &ElectronBeam::set_n_bunches);

    py::class_<UniformCylinder, ElectronBeam>(m, "UniformCylinder")
        .def(py::init<double, double, double>(),
	    py::arg("current"),
	    py::arg("radius"),
	    py::arg("neutralisation") = 2
	)
        .def("current", &UniformCylinder::current)
        .def("radius", &UniformCylinder::radius)
        .def("shape", &UniformCylinder::shape)
        .def("length", &UniformCylinder::length);

    py::class_<UniformHollow, ElectronBeam>(m, "UniformHollow")
        .def(py::init<double, double, double>())
        .def("current", &UniformHollow::current)
        .def("shape", &UniformHollow::shape)
        .def("length", &UniformHollow::length)
        .def("out_radius", &UniformHollow::out_radius)
        .def("in_radius", &UniformHollow::in_radius);

    py::class_<UniformHollowBunch, ElectronBeam>(m, "UniformHollowBunch")
        .def(py::init<double, double, double, double>())
        .def("current", &UniformHollowBunch::current)
        .def("shape", &UniformHollowBunch::shape)
        .def("length", &UniformHollowBunch::length)
        .def("out_radius", &UniformHollowBunch::out_radius)
        .def("in_radius", &UniformHollowBunch::in_radius);

    py::class_<UniformBunch, ElectronBeam>(m, "UniformBunch")
        .def(py::init<double, double, double>())
        .def("current", &UniformBunch::current)
        .def("shape", &UniformBunch::shape)
        .def("length", &UniformBunch::length)
        .def("radius", &UniformBunch::radius);

    py::class_<EllipticUniformBunch, ElectronBeam>(m, "EllipticUniformBunch")
        .def(py::init<double, double, double, double>())
        .def("shape", &EllipticUniformBunch::shape)
        .def("length", &EllipticUniformBunch::length);

    py::class_<GaussianBunch, ElectronBeam>(m, "GaussianBunch")
        .def(py::init<double, double, double, double>())
        .def("shape", &GaussianBunch::shape)
        .def("length", &GaussianBunch::length)
        .def("set_angles", &GaussianBunch::set_angles);

    py::class_<ParticleBunch, ElectronBeam>(m, "ParticleBunch")
        .def(py::init<double, std::string, double>())
        .def(py::init<double, std::string>())
        .def("shape", &ParticleBunch::shape)
        .def("length", &ParticleBunch::length)
        .def("bunched", &ParticleBunch::bunched)
        .def("corr", &ParticleBunch::corr)
        .def("set_corr", &ParticleBunch::set_corr)
        .def("set_buffer", &ParticleBunch::set_buffer)
        .def("set_s", &ParticleBunch::set_s)
        .def("set_binary", &ParticleBunch::set_binary)
        .def("set_skip", &ParticleBunch::set_skip)
        .def("load_particle", py::overload_cast<long int>(&ParticleBunch::load_particle))
        .def("load_particle", py::overload_cast<>(&ParticleBunch::load_particle));

    py::enum_<ElectronBeam::Shape>(m, "ElectronBeamShape", py::arithmetic())
        .value("UNIFORM_CYLINDER", ElectronBeam::Shape::UNIFORM_CYLINDER)
        .value("GAUSSIAN_BUNCH", ElectronBeam::Shape::GAUSSIAN_BUNCH)
        .value("UNIFORM_BUNCH", ElectronBeam::Shape::UNIFORM_BUNCH)
        .value("GAUSSIAN_CYLINDER", ElectronBeam::Shape::GAUSSIAN_CYLINDER)
        .value("ELLIPTIC_UNIFORM_BUNCH", ElectronBeam::Shape::ELLIPTIC_UNIFORM_BUNCH)
        .value("UNIFORM_HOLLO", ElectronBeam::Shape::UNIFORM_HOLLOW)
        .value("UNIFORM_HOLLOW_BUNCH", ElectronBeam::Shape::UNIFORM_HOLLOW_BUNCH)
        .value("PARTICLE_BUNCH", ElectronBeam::Shape::PARTICLE_BUNCH);

    py::enum_<ElectronBeam::Velocity>(m, "ElectronBeamVelocity", py::arithmetic())
        .value("CONST", ElectronBeam::Velocity::CONST)
        .value("USER_DEFINE", ElectronBeam::Velocity::USER_DEFINE)
        .value("SPACE_CHARGE", ElectronBeam::Velocity::SPACE_CHARGE)
        .value("VARY", ElectronBeam::Velocity::VARY)
        .value("VARY_X", ElectronBeam::Velocity::VARY_X)
        .value("VARY_Y", ElectronBeam::Velocity::VARY_Y)
        .value("VARY_Z", ElectronBeam::Velocity::VARY_Z);

    py::enum_<ElectronBeam::Temperature>(m, "ElectronBeamTemperature", py::arithmetic())
        .value("CONST", ElectronBeam::Temperature::CONST)
        .value("USER_DEFINE", ElectronBeam::Temperature::USER_DEFINE)
        .value("SPACE_CHARGE", ElectronBeam::Temperature::SPACE_CHARGE)
        .value("VARY", ElectronBeam::Temperature::VARY)
        .value("VARY_X", ElectronBeam::Temperature::VARY_X)
        .value("VARY_Y", ElectronBeam::Temperature::VARY_Y)
        .value("VARY_Z", ElectronBeam::Temperature::VARY_Z);

//    py::enum_<EdgeEffect>(m, "EdgeEffect", py::arithmetic())
//        .value("Rising", EdgeEffect::Rising)
//        .value("Falling", EdgeEffect::Falling);

}
