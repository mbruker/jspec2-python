#include <chrono>
#include <cmath>

#include "jspec2/electron_beam.h"
#include "jspec2/particle_model.h"
#include "jspec2/constants.h"
#include "jspec2/ecooling.h"
#include "jspec2/functions.h"
#include "jspec2/ion_beam.h"
#include "jspec2/cooler.h"
#include "jspec2/ring.h"

void ParticleModel::update_ibeam(IonBeam& ionBeam, ElectronBeam& ebeam, double dt) {
//    vector<double>& dp_p = ion_sample.cdnt_dp_p();
    if(ecool_solver) {
        double freq = k_c*ionBeam.beta()/ring.circ()*cooler.section_number();
        ecool_solver->adjust_rate(ionBeam, ebeam, {&freq});
        apply_cooling_kick(freq, ionBeam, dt);
    }
    if(ibs_solver) {
        apply_ibs_kick(ionBeam, dt);
    }

//    if(edge_effect) {
//        apply_edge_kick(cooler, ebeam, ionBeam);
//    }

    move_particles(ionBeam);
    update_beam_parameters(ionBeam);

    if(fixed_bunch_length && ionBeam.bunched()) {
        ring.update_bet_s();
        ring.update_rf_voltage();
    }
}

void ParticleModel::apply_cooling_kick(double freq, IonBeam& ionBeam, double dt) {
    vector<double>& xp = ionBeam.cdnt_xp();
    vector<double>& yp = ionBeam.cdnt_yp();
    vector<double>& dp_p = ionBeam.cdnt_dp_p();
    const vector<double>& force_x = ecool_solver->get_force_x();
    const vector<double>& force_y = ecool_solver->get_force_y();
    const vector<double>& force_z = ecool_solver->get_force_z();
    const double p0 = ionBeam.p0_SI();
    const double t_cooler = ecool_solver->t_cooler();
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<ionBeam.n_sample(); ++i) {
        xp[i] = !iszero(xp[i])?xp[i]*exp(force_x[i]*t_cooler*dt*freq/(xp[i]*p0)):xp[i];
        yp[i] = !iszero(yp[i])?yp[i]*exp(force_y[i]*t_cooler*dt*freq/(yp[i]*p0)):yp[i];
        const double dp = force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0);
        if(!iszero(dp_p[i],1e-7)) {
            dp_p[i] = dp>0.15?dp_p[i]*(1+dp):dp_p[i]*exp(dp);
        }
//        dp_p[i] = !iszero(dp_p[i],1e-7)?dp_p[i]*exp(force_z[i]*t_cooler*dt*freq/(dp_p[i]*p0)):dp_p[i];
    }
}

void ParticleModel::apply_ibs_kick(IonBeam& ionBeam, double dt) {
    auto twiss = ionBeam.get_twiss();
    assert(twiss.bet_x>0&& twiss.bet_y>0
           &&"TWISS parameters for the reference point not defined! Define twiss_ref.");
    ibs_kick(ionBeam.n_sample(), state.rx_ibs, twiss.bet_x, ionBeam.rms_emit_x(), ionBeam.cdnt_xp(), dt);
    ibs_kick(ionBeam.n_sample(), state.ry_ibs, twiss.bet_y, ionBeam.rms_emit_y(), ionBeam.cdnt_yp(), dt);
    ibs_kick(ionBeam.n_sample(), state.rs_ibs, ionBeam.bunched() ? 1.0 : 2.0, ionBeam.rms_dp_p()*ionBeam.rms_dp_p(), ionBeam.cdnt_dp_p(), dt);
}

void ParticleModel::ibs_kick(int n_sample, double rate, double twiss, double emit, vector<double>& v, double dt) {
    if(rate>0) {
        double theta = sqrt(2*rate*dt*emit/twiss);
        vector<double> rdn(n_sample);
        gaussian_random(n_sample, rdn, 1, 0);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) v[i] += theta*rdn[i];
    }
    else {
        double k = exp(rate*dt);
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) v[i] *= k;
    }
}

void ParticleModel::move_particles(IonBeam& ionBeam) {
    //New betatron oscillation coordinates
    auto twiss = ionBeam.get_twiss();
    vector<double>& x_bet = ionBeam.cdnt_x_bet();
    vector<double>& xp_bet = ionBeam.cdnt_xp_bet();
    vector<double>& y_bet = ionBeam.cdnt_y_bet();
    vector<double>& yp_bet = ionBeam.cdnt_yp_bet();
    vector<double>& dp_p = ionBeam.cdnt_dp_p();

    int n_sample = ionBeam.n_sample();
    ionBeam.adjust_disp_inv();

     //random phase advance
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;

    double gamma_x = (1+alf_x*alf_x)/beta_x;
    double gamma_y = (1+alf_y*alf_y)/beta_y;
    
    vector<double> rdn(n_sample);
    uniform_random(n_sample, rdn, -1, 1);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        double I = beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
        double phi = k_pi*rdn[i];
        x_bet[i] = sqrt(I*beta_x)*sin(phi);
        xp_bet[i] = sqrt(I/beta_x)*(cos(phi)-alf_x*sin(phi));
    }
    uniform_random(n_sample, rdn, -1, 1);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        double I = beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
        double phi = k_pi*rdn[i];
        y_bet[i] = sqrt(I*beta_y)*sin(phi);
        yp_bet[i] = sqrt(I/beta_y)*(cos(phi)-alf_y*sin(phi));
    }

    if(ionBeam.bunched()){
        uniform_random(n_sample, rdn, -1, 1);
        double beta_s = ring.beta_s();
        if(fixed_bunch_length) beta_s =  ionBeam.rms_sigma_s()/rms(n_sample, dp_p);
        double beta_s2_inv = 1/(beta_s*beta_s);
        vector<double>& ds = ionBeam.cdnt_ds();
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i){
            double I = ds[i]*ds[i]*beta_s2_inv+dp_p[i]*dp_p[i];
            I = sqrt(I);
            double phi = k_pi*rdn[i];
            dp_p[i] = I*sin(phi);
            ds[i] = I*beta_s*cos(phi);
        }
        if(fixed_bunch_length) {
            gaussian_random_adjust(n_sample, ds, ionBeam.rms_sigma_s());
        }
    }
    ionBeam.adjust_disp();

}

// TODO This function may need to be part of IonBeam rather than ParticleModel
void ParticleModel::update_beam_parameters(IonBeam& ionBeam) {
    const auto [emit_x, emit_y, emit_z] = ionBeam.emit();
    ionBeam.set_emit_x(emit_x);
    ionBeam.set_emit_y(emit_y);

    if(ionBeam.bunched()) {
        if(fixed_bunch_length) {
            ionBeam.set_dp_p(rms(ionBeam.n_sample(), ionBeam.cdnt_dp_p()));
            ionBeam.update_bet_s();
        }
        else {
            ionBeam.set_sigma_s(rms(ionBeam.n_sample(), ionBeam.cdnt_ds()));
            ionBeam.set_dp_p(rms(ionBeam.n_sample(), ionBeam.cdnt_dp_p()));
        }
    }
    else {
        ionBeam.set_dp_p(sqrt(emit_z));
    }
}

void ParticleModel::save_ions(double t, IonBeam& ion_sample)
{
    if (datasink)
        datasink->output_ion_phasespace(t, ion_sample);
}
