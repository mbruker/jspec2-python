#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "jspec2/datasink.h"
#include "jspec2/ion_beam.h"

namespace py=pybind11;
using namespace pybind11::literals;

struct PyPhasespace {
    double x, xp, y, yp, ds, dp_p;
    double x_bet, xp_bet, y_bet, yp_bet;
};
struct PyEcoolForce {
    double n_e, v_tr, v_l, force_x, force_z;
};
class PySink : public DataSink {
public:
    
    py::object queue;
    py::function put;
    PySink(py::object _queue);

    virtual void output_simulator_state(const Simulator::State &state) override;
    virtual void output_ion_phasespace(double t, const IonBeam &ionBeam) override;
    virtual void output_ecool_force(const vector<double> &n_e, const vector<double> &v_tr, const vector<double> &v_l, const vector<double> &force_x, const vector<double> &force_z) override;
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
    put(py::dict("type"_a="state",
                 "data"_a=state_array));
}

void PySink::output_ion_phasespace(double t, const IonBeam &ionBeam)
{
    py::array_t<PyPhasespace> ion_array(ionBeam.cdnt_x().size());
    PyPhasespace *data = static_cast<PyPhasespace *>(ion_array.request().ptr);
    for (int i = 0; i < ionBeam.cdnt_x().size(); ++i) {
        data[i].x = ionBeam.cdnt_x()[i];
        data[i].xp = ionBeam.cdnt_xp()[i];
        data[i].x_bet = ionBeam.cdnt_x_bet()[i];
        data[i].xp_bet = ionBeam.cdnt_xp_bet()[i];
        data[i].y = ionBeam.cdnt_y()[i];
        data[i].yp = ionBeam.cdnt_yp()[i];
        data[i].y_bet = ionBeam.cdnt_y_bet()[i];
        data[i].yp_bet = ionBeam.cdnt_yp_bet()[i];
        data[i].ds = ionBeam.cdnt_ds()[i];
        data[i].dp_p = ionBeam.cdnt_dp_p()[i];
    }
    put(py::dict("type"_a="phasespace",
                 "t"_a=t,
                 "data"_a=ion_array));
}

void PySink::output_ecool_force(const vector<double> &n_e, const vector<double> &v_tr, const vector<double> &v_l, const vector<double> &force_x, const vector<double> &force_z)
{
    py::array_t<PyEcoolForce> force_array(n_e.size());
    PyEcoolForce *data = static_cast<PyEcoolForce *>(force_array.request().ptr);
    for (int i = 0; i < n_e.size(); ++i) {
        data[i].n_e = n_e[i];
        data[i].v_tr = v_tr[i];
        data[i].v_l = v_l[i];
        data[i].force_x = force_x[i];
        data[i].force_z = force_z[i];
    }
    put(py::dict("type"_a="ecool_force",
                 "data"_a=force_array));
}

void init_datasink(py::module &m)
{
    py::class_<DataSink>(m, "DataSink");
    py::class_<PySink, DataSink>(m, "PySink")
        .def(py::init<py::object>());
        
    PYBIND11_NUMPY_DTYPE(Simulator::State, t, ex, ey, dp_p, sigma_s, rx, ry, rs, rx_ibs, ry_ibs, rs_ibs, rx_ecool, ry_ecool, rs_ecool, rf_voltage, luminosity);
    PYBIND11_NUMPY_DTYPE(PyPhasespace, x, xp, x_bet, xp_bet, y, yp, y_bet, yp_bet, ds, dp_p);
    PYBIND11_NUMPY_DTYPE(PyEcoolForce, n_e, v_tr, v_l, force_x, force_z);
}
