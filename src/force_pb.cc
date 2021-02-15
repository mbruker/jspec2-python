#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

#include "jspec2/force.h"

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;

void init_force(py::module& m) {
    py::class_<FrictionForceSolver>(m, "FrictionForceSolver");
    py::class_<ForceNonMag>(m, "ForceNonMag");
    py::class_<ForcePark, FrictionForceSolver>(m, "ForcePark")
        .def(py::init<>())
        .def("set_time_cooler", &ForcePark::set_time_cooler)
        .def("set_mag_field", &ForcePark::set_mag_field)
        .def("t_cooler", &ForcePark::t_cooler)
        .def("set_t_eff", &ForcePark::set_t_eff)
        .def("set_v_eff", &ForcePark::set_v_eff)
        .def("friction_force", &ForcePark::friction_force);

    py::class_<ForceNonMagDerbenev, ForceNonMag>(m, "ForceNonMagDerbenev")
        .def("set_time_cooler", &ForceNonMagDerbenev::set_time_cooler)
        .def("set_mag_field", &ForceNonMagDerbenev::set_mag_field)
        .def("t_cooler", &ForceNonMagDerbenev::t_cooler)
        .def("set_smooth_rho_max",&ForceNonMagDerbenev::set_smooth_rho_max)
        .def("friction_force",&ForceNonMagDerbenev::friction_force);

    py::class_<ForceNonMagMeshkov, ForceNonMag>(m, "ForceNonMagMeshkov")
        .def("set_time_cooler", &ForceNonMagMeshkov::set_time_cooler)
        .def("set_mag_field", &ForceNonMagMeshkov::set_mag_field)
        .def("t_cooler", &ForceNonMagMeshkov::t_cooler)
        .def("set_smooth_rho_max",&ForceNonMagMeshkov::set_smooth_rho_max)
        .def("friction_force",&ForceNonMagMeshkov::friction_force);

    py::class_<ForceNonMagNumeric1D, ForceNonMag>(m, "ForceNonMagNumeric1D")
        .def(py::init<int>())
        .def("set_time_cooler", &ForceNonMagNumeric1D::set_time_cooler)
        .def("set_mag_field", &ForceNonMagNumeric1D::set_mag_field)
        .def("t_cooler", &ForceNonMagNumeric1D::t_cooler)
        .def("set_smooth_rho_max",&ForceNonMagNumeric1D::set_smooth_rho_max)
        .def("friction_force",&ForceNonMagNumeric1D::friction_force)
        .def("set_espabs", &ForceNonMagNumeric1D::set_espabs)
        .def("set_esprel", &ForceNonMagNumeric1D::set_esprel);

    py::class_<ForceNonMagNumeric3D, ForceNonMag>(m, "ForceNonMagNumeric3D")
        .def(py::init<int>())
        .def("set_time_cooler", &ForceNonMagNumeric3D::set_time_cooler)
        .def("set_mag_field", &ForceNonMagNumeric3D::set_mag_field)
        .def("t_cooler", &ForceNonMagNumeric3D::t_cooler)
        .def("set_smooth_rho_max",&ForceNonMagNumeric3D::set_smooth_rho_max)
        .def("friction_force",&ForceNonMagNumeric3D::friction_force)
        .def("set_espabs", &ForceNonMagNumeric3D::set_espabs)
        .def("set_esprel", &ForceNonMagNumeric3D::set_esprel)
        .def("set_gsl", &ForceNonMagNumeric3D::set_gsl)
        .def("set_mean_rho_min", &ForceNonMagNumeric3D::set_mean_rho_min)
        .def("set_grid", &ForceNonMagNumeric3D::set_grid);

    py::class_<ForceMeshkov, FrictionForceSolver>(m, "ForceMeshkov")
        .def("set_time_cooler", &ForceMeshkov::set_time_cooler)
        .def("set_mag_field", &ForceMeshkov::set_mag_field)
        .def("t_cooler", &ForceMeshkov::t_cooler)
        .def("set_smooth_factor", &ForceMeshkov::set_smooth_factor)
        .def("friction_force",&ForceMeshkov::friction_force);

    py::class_<ForceDSM, FrictionForceSolver>(m, "ForceDSM")
        .def("set_time_cooler", &ForceDSM::set_time_cooler)
        .def("set_mag_field", &ForceDSM::set_mag_field)
        .def("t_cooler", &ForceDSM::t_cooler)
        .def("set_smooth_factor", &ForceDSM::set_smooth_factor)
        .def("set_steps", &ForceDSM::set_steps)
        .def("set_grid", &ForceDSM::set_grid)
        .def("set_mag_only", &ForceDSM::set_mag_only)
        .def("friction_force",&ForceDSM::friction_force);

}
