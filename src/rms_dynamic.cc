#include <chrono>
#include <cmath>

#include "jspec2/constants.h"
#include "jspec2/functions.h"
#include "jspec2/ring.h"
#include "jspec2/ion_beam.h"
#include "jspec2/rms_dynamic.h"

void RMSModel::update_ibeam(IonBeam& ionBeam, ElectronBeam& /*ebeam*/, double dt)
{
    const double emit_nx = ionBeam.rms_emit_nx() * exp(state.rx*dt);
    const double emit_ny = ionBeam.rms_emit_ny() * exp(state.ry*dt);
    // TODO Is there a reason why we compute the sqrt of a square?
    // ionBeam.rms_dp_p() can't be negative, can it?
    const double dp = sqrt(ionBeam.rms_dp_p() * ionBeam.rms_dp_p() * exp(state.rs*dt));

    ionBeam.set_emit_nx(emit_nx);
    ionBeam.set_emit_ny(emit_ny);
    ionBeam.set_dp_p(dp);

    if(ionBeam.bunched()) {
        if(fixed_bunch_length) {
            ring.update_rf_voltage();
            ionBeam.update_bet_s();
        }
        else {
            ionBeam.set_sigma_s(ionBeam.beta_s()*dp);
        }
    }

    if (ecool_solver)
        ionBeam.create_samples();
}

void RMSModel::adjust_rf_voltage()
{
    if (ring.gamma_tr() > 0)
        ring.update_rf_voltage();
}
