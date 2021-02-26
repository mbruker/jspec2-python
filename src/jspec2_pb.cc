#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <functional>
#include <tuple>
#include <vector>

namespace py=pybind11;
using namespace pybind11::literals;
using std::vector;


void init_electron_beam(py::module &);
void init_ring(py::module &);
void init_ibs(py::module &);
void init_cooler(py::module &);
void init_ecooling(py::module &);
void init_force(py::module &);
void init_ions(py::module &);
void init_luminosity(py::module &);
void init_simulators(py::module &);
void init_datasink(py::module &);

PYBIND11_MODULE(jspec, m) {
    m.doc() = "JSPEC lib";
    init_electron_beam(m);
    init_ring(m);
    init_ibs(m);
    init_cooler(m);
    init_ecooling(m);
    init_force(m);
    init_ions(m);
    init_luminosity(m);
    init_simulators(m);
    init_datasink(m);
}
