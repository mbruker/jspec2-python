#include <cmath>
#ifdef _OPENMP
    #include <omp.h>
#endif // _OPENMP

#include "jspec2/ecooling.h"
#include "jspec2/electron_beam.h"
#include "jspec2/constants.h"
#include "jspec2/cooler.h"
#include "jspec2/datasink.h"
#include "jspec2/force.h"
#include "jspec2/functions.h"
#include "jspec2/ring.h"
#include "jspec2/ion_beam.h"

void ECoolRate::electron_density(const IonBeam& ion, ElectronBeam &ebeam) {
    const int n_sample = ion.n_sample();
    const vector<double>& x = ion.cdnt_x();
    const vector<double>& y = ion.cdnt_y();
    const vector<double>& ds = ion.cdnt_ds();
    if(ebeam.p_shift()) {
        if (ebeam.multi_bunches())
            ebeam.multi_density(x, y, ds, ne, n_sample, ion.center_x(), ion.center_y(), ion.center_z());
        else
            ebeam.density(x, y, ds, ne, n_sample, ion.center_x(), ion.center_y(), ion.center_z());
    }
    else {
        if(ebeam.multi_bunches()) ebeam.multi_density(x, y, ds, ne, n_sample);
        else ebeam.density(x, y, ds, ne, n_sample);
    }
}

void ECoolRate::space_to_dynamic(const IonBeam &ion)
{
    const double v = ion.beta()*k_c;
    const vector<double>& xp = ion.cdnt_xp();
    const vector<double>& yp = ion.cdnt_yp();
    const vector<double>&dp_p = ion.cdnt_dp_p();
    const auto n_sample = xp.size();
    for(int i=0; i<n_sample; ++i) {
//        v_long[i] = dp_p[i]*v/(ion.gamma()*ion.gamma());  //Convert from dp/p to dv/v
//        v_long[i] /= (1-(v_long[i]+v)*ion.beta()/k_c);    //Convert to beam frame, when v_long<<v, canceled with the above line.
        v_long[i] = dp_p[i]*v;
        v_tr[i] = sqrt(xp[i]*xp[i]+yp[i]*yp[i])*v;
    }
}

void ECoolRate::init_scratch(int n_sample)
{
    if (n_sample != x.size()) {
        ne.resize(n_sample);
        xp_bet.resize(n_sample);
        yp_bet.resize(n_sample);
        x_bet.resize(n_sample);
        y_bet.resize(n_sample);
        xp.resize(n_sample);
        yp.resize(n_sample);
        x.resize(n_sample);
        y.resize(n_sample);
        dp_p.resize(n_sample);
        v_tr.resize(n_sample);
        v_long.resize(n_sample);
        force_x.resize(n_sample);
        force_y.resize(n_sample);
        force_z.resize(n_sample);
    }
}

void ECoolRate::beam_frame(double gamma_e)
{
    const double gamma_e_inv = 1/gamma_e;
    const auto n_sample = v_tr.size();
    for(int i=0; i<n_sample; ++i){
        v_tr[i] *= gamma_e;
        ne[i] *= gamma_e_inv;
    }
    t_cooler_ /= gamma_e;
}

void ECoolRate::lab_frame(double gamma_e)
{
    const double gamma_e_inv = 1/gamma_e;
    const auto n_sample = v_tr.size();
    for(int i=0; i<n_sample; ++i){
            force_x[i] *= gamma_e_inv;
            v_tr[i] *= gamma_e_inv;
    }
    t_cooler_ *= gamma_e;
}

//Calculate friction force
void ECoolRate::force(const IonBeam &ion, const ElectronBeam &ebeam, const Cooler &cooler)
{
    const auto n_sample = ion.n_sample();
    
    //set parameters for friction force calculation
    force_solver->set_mag_field(cooler.magnetic_field());
    force_solver->set_time_cooler(t_cooler_);
    const double d_beta = ebeam.beta() - ion.beta();
    double cv_l = 0;
    if(!iszero(d_beta, 1e-6)) cv_l = d_beta*k_c;
    if(!iszero(ebeam.cv_l(), 1e-6)) cv_l += ebeam.cv_l();

    if(!iszero(cv_l, 1e-6)) {
        for(auto& v: v_long) {
            v -= cv_l;
        }
    }
    force_solver->friction_force(ion.charge_number(), n_sample, v_tr, v_long, ne, ebeam, force_x, force_z);

    if(longitudinal_force_solver) {
        longitudinal_force_solver->set_mag_field(cooler.magnetic_field());
        longitudinal_force_solver->set_time_cooler(t_cooler_);
        longitudinal_force_solver->friction_force(ion.charge_number(), n_sample, v_tr, v_long, ne, ebeam, force_y, force_z); //force_y will be ignored in the following.
    }

    if(datasink)
        datasink->output_ecool_force(ne, v_tr, v_long, force_x, force_z);
}

void ECoolRate::bunched_to_coasting(IonBeam& ion, ElectronBeam &ebeam, const Cooler &cooler)
{
    const int n_sample = ion.n_sample();
    int count = 1;
    vector<double> force_tr_rcd(n_sample);
    vector<double> force_long_rcd(n_sample);

    const double cz_rcd = ion.center_z();

    ebeam.set_p_shift(true);
    const double length = ebeam.length();
    const double step = length/n_long_sample_;
    const double gamma_e_inv = 1/ebeam.gamma();
    for(double cz = cz_rcd-0.5*length; cz <= cz_rcd+0.5*length; cz += step) {
        ion.set_center_z(cz);
        electron_density(ion, ebeam);
        for(int i=0; i<n_sample; ++i) ne[i] *= gamma_e_inv;
        force(ion, ebeam, cooler);
        for(int i=0; i<n_sample; ++i) {
            force_tr_rcd[i] += force_x[i];
            force_long_rcd[i] += force_z[i];
        }
        ++count;
    }

    ion.set_center_z(cz_rcd);
    ebeam.set_p_shift(false);
    const double count_inv = 1.0/static_cast<double>(count);
    for(int i=0; i<n_sample; ++i) {
        force_x[i] = force_tr_rcd[i]*count_inv;
        force_z[i] = force_long_rcd[i]*count_inv;
    }
}

//Distribute to x and y direction
void ECoolRate::force_distribute(const IonBeam &ion)
{
    const double v0 = ion.beta()*k_c;
    const vector<double>& xp = ion.cdnt_xp();
    const vector<double>& yp = ion.cdnt_yp();
    const auto n_sample = xp.size();
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        force_y[i] = yp[i]!=0?force_x[i]*yp[i]*v0/v_tr[i]:0;
        force_x[i] = xp[i]!=0?force_x[i]*xp[i]*v0/v_tr[i]:0;
    }
}

void ECoolRate::apply_kick(const IonBeam& ion)
{
    const double p0 = ion.p0_SI();
    const vector<double>& ixp = ion.cdnt_xp();
    const vector<double>& iyp = ion.cdnt_yp();
    const vector<double>& idp_p = ion.cdnt_dp_p();
    const auto n_sample = ixp.size();
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n_sample; ++i){
        xp[i] = !iszero(ixp[i])?ixp[i]*exp(force_x[i]*t_cooler_/(p0*ixp[i])):ixp[i];
        yp[i] = !iszero(iyp[i])?iyp[i]*exp(force_y[i]*t_cooler_/(p0*iyp[i])):iyp[i];
//        dp_p[i] = !iszero(idp_p[i])?idp_p[i]*exp(force_z[i]*t_cooler_/(p0*idp_p[i])):idp_p[i];
        double dp = force_z[i]*t_cooler_/(idp_p[i]*p0);
        if(!iszero(idp_p[i],1e-7)) {
            dp_p[i] = dp>0.15?idp_p[i]*(1+dp):idp_p[i]*exp(dp);
        }
        else {
            dp_p[i] = idp_p[i];
        }
    }
}

void ECoolRate::adjust_rate(const IonBeam &ion, const ElectronBeam &ebeam, initializer_list<double*> func) {
    if(ebeam.bunched()&&(!ion.bunched())) {
        double sample_length = ebeam.length();
        if(sample_length<0) perror("electron bunch length must be positive!");
        if(bunch_separate_>=sample_length) {
            double duty_factor = sample_length/bunch_separate_;
            for(auto& f: func) *f *= duty_factor;
        }
        else {
            perror("Electron bunch length is larger than the distance between electron bunches");
        }
    }
}

rate3d ECoolRate::ecool_rate(IonBeam &ion, const Cooler &cooler, ElectronBeam &ebeam, const Ring &ring)
{
    const auto n_sample = ion.n_sample();
    init_scratch(n_sample);

    electron_density(ion, ebeam);
    space_to_dynamic(ion);

    //Time through the cooler
    t_cooler_ = cooler.length()/(ion.beta()*k_c);
    //Transfer into e- beam frame
    beam_frame(ebeam.gamma());
    //Calculate friction force
    force(ion, ebeam, cooler);
//    //Restore the longitudinal velocity if it has been changed due to electron velocity gradient
//    restore_velocity(ebeam);

    //Special treatment for bunched electron beam to cool coasting ion beam
    if(!ion.bunched()&&ebeam.bunched()) {
        bunched_to_coasting(ion, ebeam, cooler);
        force(ion, ebeam, cooler);
    }

    //Transfer back to lab frame
    lab_frame(ebeam.gamma());

    //Distribute the friction force into x,y direction.
    force_distribute(ion);

    //Original emittance
    const auto [emit_x0, emit_y0, emit_z0] = ion.emit();

    //Apply kick
    apply_kick(ion);

    //New emittance
    auto t = ion.get_twiss();

    adjust_disp_inv(t.disp_x, x_bet, dp_p, ion.cdnt_x(), n_sample);
    adjust_disp_inv(t.disp_dx, xp_bet, dp_p, xp, n_sample);
    adjust_disp_inv(t.disp_y, y_bet, dp_p, ion.cdnt_y(), n_sample);
    adjust_disp_inv(t.disp_dy, yp_bet, dp_p, yp, n_sample);

    const auto [emit_x, emit_y, emit_z] = ion.emit(x_bet, xp_bet, y_bet, yp_bet, dp_p, ion.cdnt_ds());

    double rate_x = emit_x/emit_x0-1;
    double rate_y = emit_y/emit_y0-1;
    double rate_s = emit_z/emit_z0-1;
    double freq = k_c*ion.beta()/ring.circ()*cooler.section_number();
    rate_x *= freq;
    rate_y *= freq;
    rate_s *= freq;

    adjust_rate(ion, ebeam, {&rate_x, &rate_y, &rate_s});
    return std::make_tuple(rate_x, rate_y, rate_s);
}


void ForceCurve::output_force(const IonBeam &ionBeam, const Cooler &cooler, ElectronBeam &ebeam)
{
    // TODO Something evil is happening in this function. A second ion beam?
    
    if(iszero(dp_p, 1e-14)) n_l = 0;
    if(iszero(angle, 1e-14)) n_tr = 0;
    const int n_sample = (n_tr+1) * (2*n_l+1);
    init_scratch(n_sample);
    IonBeam_MonteCarlo ion(cooler.twiss(), ionBeam.charge_number(), ionBeam.mass(), ionBeam.kinetic_energy(), ionBeam.rms_emit_nx(), ionBeam.rms_emit_ny(), ionBeam.rms_dp_p(), ionBeam.rms_sigma_s(), ionBeam.particle_number(), n_sample);
    if(density_e>0) ne.assign(n_sample, density_e*ebeam.gamma());
    else electron_density(ion,ebeam);

    // space to dynamic
    const double v = ionBeam.beta()*k_c;
    if (dp_p<0) dp_p *= -1;
    if (angle<0) angle *= -1;
    const double da = angle/n_tr;
    const double dp = v*dp_p/n_l;
    double dpp = dp_p*v;
    for(int i=0; i<2*n_l+1; ++i) {
        double a = angle;
        for(int j=0; j<n_tr+1; ++j) {
            v_long[i*(n_tr+1)+j] = dpp;
            v_tr[i*(n_tr+1)+j] = v*sin(a);
            a -= da;
        }
        dpp -= dp;
    }

    //Time through the cooler
    t_cooler_ = cooler.length()/v;
    // Transfer into e- beam frame
    beam_frame(ebeam.gamma());
    //Calculate friction force
    force(ionBeam, ebeam, cooler);
}

