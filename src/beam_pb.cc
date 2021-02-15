#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/beam.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_beam(py::module &m) {
    py::class_<Beam>(m, "Beam")
        .def("set_emit_nx", (void (Beam::*)(double)) &Beam::set_emit_nx, "x"_a)
        .def("set_emit_ny", (void (Beam::*)(double)) &Beam::set_emit_ny, "x"_a)
        .def("set_emit_x", (void (Beam::*)(double)) &Beam::set_emit_x, "x"_a)
        .def("set_emit_y", (void (Beam::*)(double)) &Beam::set_emit_y, "x"_a)
        .def("set_dp_p", (void (Beam::*)(double)) &Beam::set_dp_p, "x"_a)
        .def("set_center", (void (Beam::*)(double, double, double)) &Beam::set_center, "cx"_a, "cy"_a, "cz"_a)
        .def("set_center", (void (Beam::*)(int, double)) &Beam::set_center, "i"_a, "x"_a)
        .def("charge_number", &Beam::charge_number)
        .def("mass", &Beam::mass)
        .def("kinetic_energy", &Beam::kinetic_energy)
        .def("beta", &Beam::beta)
        .def("gamma", &Beam::gamma)
        .def("emit_nx", &Beam::emit_nx)
        .def("emit_ny", &Beam::emit_ny)
        .def("emit_x", &Beam::emit_x)
        .def("emit_y", &Beam::emit_y)
        .def("dp_p", &Beam::dp_p)
        .def("energy_spread", &Beam::energy_spread)
        .def("sigma_s", &Beam::sigma_s)
        .def("r", &Beam::r)
        .def("particle_number", &Beam::particle_number)
//        .def("mass_number", &Beam::mass_number)
        .def("p0_SI", &Beam::p0_SI)
        .def("p0", &Beam::p0)
        .def("bunched", &Beam::bunched)
        .def("center", (int (Beam::*)(double&, double&, double&)) &Beam::center, "cx"_a, "cy"_a, "cz"_a)
        .def("center", (double (Beam::*)(int)) &Beam::center, "i"_a)
        .def(py::init<int, double, double, double, double, double, double, double>())
        .def(py::init<int, double, double, double, double, double, double>());

    py::class_<EBeam>(m, "EBeam")
        .def("velocity", &EBeam::velocity)
        .def("temperature", &EBeam::temperature)
        .def("charge_number", &EBeam::charge_number)
        .def("mass", &EBeam::mass)
        .def("mass_number", &EBeam::mass_number)
        .def("mass_SI", &EBeam::mass_SI)
        .def("kinetic_energy", &EBeam::kinetic_energy)
        .def("gamma", &EBeam::gamma)
        .def("beta", &EBeam::beta)
        .def("bunched", &EBeam::bunched)
        .def("set_p_shift", &EBeam::set_p_shift)
        .def("set_v_shift", &EBeam::set_v_shift)
        .def("p_shift", &EBeam::p_shift)
        .def("v_shift", &EBeam::v_shift)
        .def("set_cv_l", &EBeam::set_cv_l)
        .def("cv_l", &EBeam::cv_l)
        .def("shape", &EBeam::shape)
        .def("length", &EBeam::length)
        .def("neutral", &EBeam::neutral)
        .def("set_kinetic_energy", &EBeam::set_kinetic_energy)
        .def("set_gamma", &EBeam::set_gamma)
        .def("set_beta", &EBeam::set_beta)
        .def("set_center", (int (EBeam::*)(double, double, double)) &EBeam::set_center, "cx"_a, "cy"_a, "cz"_a)
        .def("set_center", (void (EBeam::*)(int, double)) &EBeam::set_center, "i"_a, "x"_a)
        .def("center", py::overload_cast<double&, double&, double&>(&EBeam::center, py::const_))
        .def("center",  py::overload_cast<int>(&EBeam::center, py::const_))
        .def("set_tpr", &EBeam::set_tpr)
        .def("set_v_rms", &EBeam::set_v_rms)
        .def("set_v_avg", &EBeam::set_v_avg)
        .def("set_neutral", &EBeam::set_neutral)
        .def("set_multi_bunches", &EBeam::set_multi_bunches)
        .def("multi_bunches", &EBeam::multi_bunches)
        .def("get_v", &EBeam::get_v)
        .def("set_n_bunches", &EBeam::set_n_bunches);

    py::class_<UniformCylinder, EBeam>(m, "UniformCylinder")
        .def(py::init<double, double, double>(),
	    py::arg("current"),
	    py::arg("radius"),
	    py::arg("neutralisation") = 2
	)
        .def("current", &UniformCylinder::current)
        .def("radius", &UniformCylinder::radius)
        .def("shape", &UniformCylinder::shape)
        .def("length", &UniformCylinder::length);

    py::class_<UniformHollow, EBeam>(m, "UniformHollow")
        .def(py::init<double, double, double>())
        .def("current", &UniformHollow::current)
        .def("shape", &UniformHollow::shape)
        .def("length", &UniformHollow::length)
        .def("out_radius", &UniformHollow::out_radius)
        .def("in_radius", &UniformHollow::in_radius);

    py::class_<UniformHollowBunch, EBeam>(m, "UniformHollowBunch")
        .def(py::init<double, double, double, double>())
        .def("current", &UniformHollowBunch::current)
        .def("shape", &UniformHollowBunch::shape)
        .def("length", &UniformHollowBunch::length)
        .def("out_radius", &UniformHollowBunch::out_radius)
        .def("in_radius", &UniformHollowBunch::in_radius);

    py::class_<UniformBunch, EBeam>(m, "UniformBunch")
        .def(py::init<double, double, double>())
        .def("current", &UniformBunch::current)
        .def("shape", &UniformBunch::shape)
        .def("length", &UniformBunch::length)
        .def("radius", &UniformBunch::radius);

    py::class_<EllipticUniformBunch, EBeam>(m, "EllipticUniformBunch")
        .def(py::init<double, double, double, double>())
        .def("shape", &EllipticUniformBunch::shape)
        .def("length", &EllipticUniformBunch::length);

    py::class_<GaussianBunch, EBeam>(m, "GaussianBunch")
        .def(py::init<double, double, double, double>())
        .def("shape", &GaussianBunch::shape)
        .def("length", &GaussianBunch::length)
        .def("set_angles", &GaussianBunch::set_angles);

    py::class_<ParticleBunch, EBeam>(m, "ParticleBunch")
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

    py::enum_<Shape>(m, "Shape", py::arithmetic())
        .value("UNIFORM_CYLINDER", Shape::UNIFORM_CYLINDER)
        .value("GAUSSIAN_BUNCH", Shape::GAUSSIAN_BUNCH)
        .value("UNIFORM_BUNCH", Shape::UNIFORM_BUNCH)
        .value("GAUSSIAN_CYLINDER", Shape::GAUSSIAN_CYLINDER)
        .value("ELLIPTIC_UNIFORM_BUNCH", Shape::ELLIPTIC_UNIFORM_BUNCH)
        .value("UNIFORM_HOLLO", Shape::UNIFORM_HOLLOW)
        .value("UNIFORM_HOLLOW_BUNCH", Shape::UNIFORM_HOLLOW_BUNCH)
        .value("PARTICLE_BUNCH", Shape::PARTICLE_BUNCH);

    py::enum_<Velocity>(m, "Velocity", py::arithmetic())
        .value("CONST", Velocity::CONST)
        .value("USER_DEFINE", Velocity::USER_DEFINE)
        .value("SPACE_CHARGE", Velocity::SPACE_CHARGE)
        .value("VARY", Velocity::VARY)
        .value("VARY_X", Velocity::VARY_X)
        .value("VARY_Y", Velocity::VARY_Y)
        .value("VARY_Z", Velocity::VARY_Z);

    py::enum_<Temperature>(m, "Temperature", py::arithmetic())
        .value("CONST", Temperature::CONST)
        .value("USER_DEFINE", Temperature::USER_DEFINE)
        .value("SPACE_CHARGE", Temperature::SPACE_CHARGE)
        .value("VARY", Temperature::VARY)
        .value("VARY_X", Temperature::VARY_X)
        .value("VARY_Y", Temperature::VARY_Y)
        .value("VARY_Z", Temperature::VARY_Z);

    py::enum_<EBeamV>(m, "EBeamV", py::arithmetic())
        .value("TPR_TR", EBeamV::TPR_TR)
        .value("TPR_L", EBeamV::TPR_L)
        .value("V_RMS_TR", EBeamV::V_RMS_TR)
        .value("V_RMS_L", EBeamV::V_RMS_L)
        .value("V_AVG_X", EBeamV::V_AVG_X)
        .value("V_AVG_Y", EBeamV::V_AVG_Y)
        .value("V_AVG_L", EBeamV::V_AVG_L);

//    py::enum_<EdgeEffect>(m, "EdgeEffect", py::arithmetic())
//        .value("Rising", EdgeEffect::Rising)
//        .value("Falling", EdgeEffect::Falling);

}
