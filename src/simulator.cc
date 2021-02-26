#include <cmath>

#include "jspec2/simulator.h"
#include "jspec2/electron_beam.h"
#include "jspec2/ring.h"
#include "jspec2/ibs.h"
#include "jspec2/ecooling.h"
#include "jspec2/functions.h"
#include "jspec2/luminosity.h"
#include "jspec2/ion_beam.h"
#include "jspec2/datasink.h"

double Simulator::calc_timestep(double time, int n_steps) const
{
    return time / n_steps;
}

void Simulator::run(IonBeam& ionBeam, ElectronBeam& ebeam,
                    double time, int n_steps, int state_output_interval, int ion_output_interval)
{
    double dt = calc_timestep(time, n_steps);
    if (!n_steps)
        n_steps = ceil(time / dt);
    
    adjust_rf_voltage();

    for(int i=0; i<=n_steps; ++i) {
        if (i % ion_output_interval == 0)
            save_ions(state.t, ionBeam);
        
        state.ex = ionBeam.rms_emit_nx();
        state.ey = ionBeam.rms_emit_ny();
        state.dp_p = ionBeam.rms_dp_p();
        if (ionBeam.bunched())
            state.sigma_s = ionBeam.rms_sigma_s();

        if(ibs_solver) {
            std::tie(state.rx_ibs, state.ry_ibs, state.rs_ibs) = ibs_solver->rate(ring.lattice(), ionBeam);
        }

        if(ecool_solver) {
            std::tie(state.rx_ecool, state.ry_ecool, state.rs_ecool) = ecool_solver->ecool_rate(ionBeam, cooler, ebeam, ring);
        }

        state.rx = state.rx_ibs + state.rx_ecool;
        state.ry = state.ry_ibs + state.ry_ecool;
        state.rs = state.rs_ibs + state.rs_ecool;
        state.rf_voltage = ring.rf_voltage();

        if (i%state_output_interval==0) {
            if (lum_solver) {
                if (lum_solver->use_ion_emit())
                    lum_solver->set_geo_emit(ionBeam.rms_emit_x(), ionBeam.rms_emit_y(), 0);
                state.luminosity = lum_solver->luminosity() * 10000.0;
            }
            if (datasink)
                datasink->output_simulator_state(state);
        }

        update_ibeam(ionBeam, ebeam, dt);

        state.t += dt;
    }
    save_ions(state.t, ionBeam);

}
