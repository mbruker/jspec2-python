#include <cmath>

#include "jspec2/simulator.h"
#include "jspec2/beam.h"
#include "jspec2/ring.h"
#include "jspec2/ibs.h"
#include "jspec2/ecooling.h"
#include "jspec2/functions.h"
#include "jspec2/luminosity.h"

double Simulator::calc_timestep(double time, int n_steps) const
{
    return time / n_steps;
}

void Simulator::run(Beam& ion, Ions& ion_sample, EBeam& ebeam,
                    double time, int n_steps, int state_output_interval, int ion_output_interval)
{
    double dt = calc_timestep(time, n_steps);
    if (!n_steps)
        n_steps = ceil(time / dt);
    
    std::cout << "Starting: dt=" << dt << "; n_steps=" << n_steps << std::endl;
    
    //Set time for new simulation or continued simulation.
//    if(reset_time) t = t0;
//    else t = uircd.t;
    
//    output_to_file();

    adjust_rf_voltage();

    State state;
    
    for(int i=0; i<=n_steps; ++i) {
//        save_ions(i, ion_sample);
        
        state.ex = ion.emit_nx();
        state.ey = ion.emit_ny();
        state.dp_p = ion.dp_p();
        if (ion.bunched())
            state.sigma_s = ion.sigma_s();

        if(ibs_solver) {
            std::tie(state.rx_ibs, state.ry_ibs, state.rs_ibs) = ibs_solver->rate(ring.lattice(), ion);
        }

        if(ecool_solver) {
            std::tie(state.rx_ecool, state.ry_ecool, state.rs_ecool) = ecool_solver->ecool_rate(*force_solver, ion, ion_sample, cooler, ebeam, ring);
        }

        state.rx = state.rx_ibs + state.rx_ecool;
        state.ry = state.ry_ibs + state.ry_ecool;
        state.rs = state.rs_ibs + state.rs_ecool;
        state.rf_voltage = ring.rf_voltage();

        if (i%state_output_interval==0) {
            if (lum_solver) {
                if (lum_solver->use_ion_emit())
                    lum_solver->set_geo_emit(ion.emit_x(), ion.emit_y(), 0);
                state.luminosity = lum_solver->luminosity() * 10000.0;
            }
            if (datasink)
                datasink->output_simulator_state(state);
        }

        update_ibeam(ion, ion_sample, ebeam, dt);

        state.t += dt;
    }
//    save_ions(n_steps, ion_sample);

}
