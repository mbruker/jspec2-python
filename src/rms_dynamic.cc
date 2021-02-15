#include <chrono>
#include <cmath>
#include "jspec2/constants.h"
#include "jspec2/functions.h"
#include "jspec2/ring.h"
#include "jspec2/ions.h"
#include "jspec2/rms_dynamic.h"

void RMSModel::update_ibeam(Beam& ion, Ions& ion_sample, EBeam& ebeam, double dt) {
    double emit_nx = ion.emit_nx();
    double emit_ny = ion.emit_ny();
    double dp = ion.dp_p();
    emit_nx *= exp(state.rx*dt);
    emit_ny *= exp(state.ry*dt);
    dp *= dp*exp(state.rs*dt);
    dp = sqrt(dp);

    ion.set_emit_nx(emit_nx);
    ion.set_emit_ny(emit_ny);
    ion.set_dp_p(dp);

    if(ion.bunched()) {
        if(fixed_bunch_length) {
            ring.update_rf_voltage();
            ring.update_bet_s();
        }
        else {
            ion.set_sigma_s(ring.beta_s()*dp);
        }
    }

    if(ecool_solver) ion_sample.create_samples(ion);
}

void RMSModel::adjust_rf_voltage()
{
    if (ring.gamma_tr() > 0)
        ring.update_rf_voltage();
}
