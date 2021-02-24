#include <assert.h>
#include <cmath>
#include <fstream>
#include <chrono>

#include "jspec2/constants.h"
#include "jspec2/ion_beam.h"
#include "jspec2/arbitrary_electron_beam.h"
#include "jspec2/functions.h"

// Changed mass to MeV/c^2 to match the legacy input file format!
IonBeam::IonBeam(const Twiss &_twiss, int charge_number, double mass, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle)
    : twiss(_twiss),
      charge_number_(charge_number),
      mass_(mass),
      kinetic_energy_(kinetic_energy),
      rms_emit_nx_(emit_nx),
      rms_emit_ny_(emit_ny),
      rms_dp_p_(dp_p),
      rms_sigma_s_(sigma_s),
      particle_number_(n_particle)
{
    gamma_ = 1+kinetic_energy_/mass_;
    beta_ = sqrt(gamma_*gamma_-1)/gamma_;
    r_ = k_ke*charge_number_*charge_number_*k_e*1e-6/mass_;
    bunched_ = (rms_sigma_s_>0)?true:false;
    rms_emit_x_ = rms_emit_nx_/(beta_*gamma_);
    rms_emit_y_ = rms_emit_ny_/(beta_*gamma_);
    rms_energy_spread_ = beta_*beta_*rms_dp_p_;
    rms_dv_v_ = rms_dp_p_/(gamma_*gamma_);
    p0_SI_ = gamma_*mass_*1e6*k_e*beta_/k_c;
}

//Calculate the transverse emittance statistically
double IonBeam_MonteCarlo::statistical_emittance_tr(const vector<double>& x, const vector<double>&xp) const
{
    const auto n = x.size();
    double x_mean = 0, xp_mean = 0, dlt2_x = 0, dlt2_xp = 0, dlt_xxp = 0;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:x_mean,xp_mean)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        x_mean += x[i];
        xp_mean += xp[i];
    }
    x_mean /= n;
    xp_mean /= n;

    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:dlt2_x,dlt2_xp,dlt_xxp)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        const double x_adj = x[i]-x_mean;
        const double xp_adj = xp[i]-xp_mean;
        dlt2_x += x_adj*x_adj;
        dlt2_xp += xp_adj*xp_adj;
        dlt_xxp += x_adj*xp_adj;
    }
    return sqrt(dlt2_x*dlt2_xp-dlt_xxp*dlt_xxp)/n;
}

//Calculate the longitudinal emittance as (dp/p)^2/n
double IonBeam_MonteCarlo::statistical_emittance_l(const vector<double>& dp_p) const
{
    const auto n = dp_p.size();
    double emit_p = 0;
    double dp_p_mean = 0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:dp_p_mean)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        dp_p_mean += dp_p[i];
    }
    dp_p_mean /= n;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:emit_p)
    #endif // _OPENMP
    for(int i=0; i<n; ++i){
        const double dp_p_adj = dp_p[i] - dp_p_mean;
        emit_p += dp_p_adj*dp_p_adj;
    }
    emit_p /= n;
    return emit_p;
}

void adjust_disp_inv(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n) {
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for(int i=0; i<n; ++i) x_bet[i] = x[i]-dx*dp_p[i];
}

void adjust_disp(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n){
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) x[i] = x_bet[i]+dx*dp_p[i];
}

array<double,3> IonBeam_MonteCarlo::emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>& ds) const
{
    double emit_x = statistical_emittance_tr(x_bet, xp_bet);
    double emit_y = statistical_emittance_tr(y_bet, yp_bet);
    double emit_s = statistical_emittance_l(dp_p);
    if(bunched_)
        emit_s += statistical_emittance_l(ds)/(beta_s_*beta_s_);
    return {emit_x, emit_y, emit_s};
}

array<double,3> IonBeam_MonteCarlo::emit() const
{
    return emit(x_bet, xp_bet, y_bet, yp_bet, dp_p, ds);
}

array<double,3> IonBeam_SingleParticle::emit(const vector<double>& x_bet,
                                          const vector<double>& xp_bet,
                                          const vector<double>& y_bet,
                                          const vector<double>& yp_bet,
                                          const vector<double>& dp_p,
                                          const vector<double>& ds) const
{
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;
    double gamma_x = (1+alf_x*alf_x)/beta_x;
    double gamma_y = (1+alf_y*alf_y)/beta_y;

    double emit_x = 0;
    double emit_y = 0;
    double emit_s = 0;
    int n_sample = n_;
    double inv_beta_s2 = 0;
    if(bunched_) inv_beta_s2 = 1/(beta_s_*beta_s_);
    for(int i=0; i<n_sample; ++i) {
        emit_x += beta_x*xp_bet[i]*xp_bet[i]+2*alf_x*x_bet[i]*xp_bet[i]+gamma_x*x_bet[i]*x_bet[i];
        emit_y += beta_y*yp_bet[i]*yp_bet[i]+2*alf_y*y_bet[i]*yp_bet[i]+gamma_y*y_bet[i]*y_bet[i];
        emit_s += dp_p[i]*dp_p[i];
        if(bunched_) emit_s += ds[i]*ds[i]*inv_beta_s2;
    }
    emit_x /= 2*n_sample;
    emit_y /= 2*n_sample;
    emit_s /= n_sample;
    
    return {emit_x, emit_y, emit_s};
}

array<double,3> IonBeam_SingleParticle::emit() const
{
    return emit(x_bet, xp_bet, y_bet, yp_bet, dp_p, ds);
}

//Generate Gaussian random number in S frame with given Twiss parameters
//First, rotate to O frame where alf = 0;
//Second, Generate x and xp with Gaussian random number in O frame
//Third, rotate back to S frame
int gaussian_bet_cod(double beta_xs, double alf_xs, double emit_x, vector<double>& x_bet, vector<double>& xp_bet, int n){

    double gamma_xs = (1+alf_xs*alf_xs)/beta_xs;
    double theta = atan(2*alf_xs/(gamma_xs-beta_xs))/2;     //rotation angle between O frame and S frame

    //Transfer matrix between O and S frames
    double matrix_os[2][2], matrix_so[2][2];
    matrix_os[0][0] = cos(theta);
    matrix_os[0][1] = -sin(theta);
    matrix_os[1][0] = sin(theta);
    matrix_os[1][1] = cos(theta);
    matrix_so[0][0] = matrix_os[0][0];
    matrix_so[0][1] = -matrix_os[0][1] ;
    matrix_so[1][0] = -matrix_os[1][0];
    matrix_so[1][1] = matrix_os[1][1];

    //Calculate beta and sigma in O frame
    double beta_xo = matrix_so[0][0]*matrix_so[0][0] * beta_xs-2*matrix_so[0][0]*matrix_so[0][1]*alf_xs+
                     matrix_so[0][1]*matrix_so[0][1]*gamma_xs;
    double sigma_xo = sqrt(emit_x*beta_xo);
    double sigma_xpo = sqrt(emit_x/beta_xo);
    //Generate x and xp in O frame
    gaussian_random(n, x_bet, sigma_xo);
    gaussian_random(n, xp_bet, sigma_xpo);
    gaussian_random_adjust(n, x_bet, sigma_xo);
    gaussian_random_adjust(n, xp_bet, sigma_xpo);

    //Rotate back to S frame
    for(int i=0; i<n;++i){
        double x = matrix_os[0][0]*x_bet[i]+matrix_os[0][1]*xp_bet[i];
        double xp = matrix_os[1][0]*x_bet[i]+matrix_os[1][1]*xp_bet[i];
        x_bet[i] = x;
        xp_bet[i] = xp;
    }
    return 0;
}

void IonBeam::adjust_disp(){
    ::adjust_disp(twiss.disp_x, x_bet, dp_p, x, n_);
    ::adjust_disp(twiss.disp_y, y_bet, dp_p, y, n_);
    ::adjust_disp(twiss.disp_dx, xp_bet, dp_p, xp, n_);
    ::adjust_disp(twiss.disp_dy, yp_bet, dp_p, yp, n_);
}

void IonBeam::adjust_disp_inv(){
    ::adjust_disp_inv(twiss.disp_x, x_bet, dp_p, x, n_);
    ::adjust_disp_inv(twiss.disp_y, y_bet, dp_p, y, n_);
    ::adjust_disp_inv(twiss.disp_dx, xp_bet, dp_p, xp, n_);
    ::adjust_disp_inv(twiss.disp_dy, yp_bet, dp_p, yp, n_);
}

IonBeam_MonteCarlo::IonBeam_MonteCarlo(const Twiss &_twiss, int _charge_number, double _mass, double _kinetic_energy, double _emit_nx, double _emit_ny, double _dp_p,
           double _sigma_s, double _n_particle, int n_sample)
    : IonBeam(_twiss, _charge_number, _mass, _kinetic_energy, _emit_nx, _emit_ny, _dp_p, _sigma_s, _n_particle)
{
    n_=n_sample;
    x_bet.resize(n_sample,0);
    y_bet.resize(n_sample,0);
    xp_bet.resize(n_sample,0);
    yp_bet.resize(n_sample,0);
    ds.resize(n_sample,0);
    dp_p.resize(n_sample,0);
    x.resize(n_sample,0);
    y.resize(n_sample,0);
    xp.resize(n_sample,0);
    yp.resize(n_sample,0);
}

IonBeam_MonteCarlo::IonBeam_MonteCarlo(const Twiss &_twiss, int _charge_number, double _mass, double _kinetic_energy, double _emit_nx, double _emit_ny, double _dp_p,
           double _sigma_s, double _n_particle, std::string filename, int n, int skip, bool binary, int n_buffer)
    : IonBeam(_twiss, _charge_number, _mass, _kinetic_energy, _emit_nx, _emit_ny, _dp_p, _sigma_s, _n_particle)
{
    auto n_loaded = load_electrons(x, xp, y, yp, ds, dp_p, filename, n, skip, binary,n_buffer);
    if (n_loaded!=n_) n_ = n_loaded;
}

void IonBeam_MonteCarlo::create_samples()
{
    gaussian_bet_cod(twiss.bet_x, twiss.alf_x, rms_emit_x_, x_bet, xp_bet, n_);
    gaussian_bet_cod(twiss.bet_y, twiss.alf_y, rms_emit_y_, y_bet, yp_bet, n_);

    gaussian_random(n_, dp_p, rms_dp_p_);
    gaussian_random_adjust(n_, dp_p, rms_dp_p_);

    //longitudinal sampling
    if(bunched()) {
        gaussian_random(n_, ds, rms_sigma_s_);
        gaussian_random_adjust(n_, ds, rms_sigma_s_);
    }

    ::adjust_disp(twiss.disp_x, x_bet, dp_p, x, n_);
    ::adjust_disp(twiss.disp_y, y_bet, dp_p, y, n_);
    ::adjust_disp(twiss.disp_dx, xp_bet, dp_p, xp, n_);
    ::adjust_disp(twiss.disp_dy, yp_bet, dp_p, yp, n_);

    if(bunched())
        update_bet_s();
}

IonBeam_SingleParticle::IonBeam_SingleParticle(const Twiss &_twiss, int _charge_number, double _mass, double _kinetic_energy, double _emit_nx, double _emit_ny, double _dp_p,
           double _sigma_s, double _n_particle, int n_tr, int n_l)
    : IonBeam(_twiss, _charge_number, _mass, _kinetic_energy, _emit_nx, _emit_ny, _dp_p, _sigma_s, _n_particle),
      n_tr_(n_tr),
      n_l_(n_l)
{
    n_=n_tr*n_tr*n_l;
    x_spl.resize(n_tr_);
    y_spl.resize(n_tr_);
    xp_spl.resize(n_tr_);
    yp_spl.resize(n_tr_);
    dp_p_spl.resize(n_l_);
    ds_spl.resize(n_l_);

    x_bet.resize(n_);
    y_bet.resize(n_);
    xp_bet.resize(n_);
    yp_bet.resize(n_);
    ds.resize(n_);
    dp_p.resize(n_);
    x.resize(n_);
    y.resize(n_);
    xp.resize(n_);
    yp.resize(n_);
};

void IonBeam_SingleParticle::single_particle_grid()
{
    const double alf_x = twiss.alf_x;
    const double alf_y = twiss.alf_y;
    double dphi = 2.0*k_pi/n_tr_;
    double phi = 0;
    for(int i=0; i<n_tr_; ++i){
        x_spl[i] = sin(phi);
        y_spl[i] = sin(phi);
        xp_spl[i] = cos(phi)-alf_x*sin(phi);
        yp_spl[i] = cos(phi)-alf_y*sin(phi);
        phi += dphi;
    }

    if(bunched()){
        phi = 0;
        dphi = 2.0*k_pi/n_l_;
        for(int i=0; i<n_l_; ++i){
            ds_spl[i] = sin(phi);
            dp_p_spl[i] = cos(phi);
            phi += dphi;
        }
    }
}

void IonBeam_SingleParticle::create_samples()
{
    // TODO eliminate redundant local variables
    double emit_x = rms_emit_x_;
    double emit_y = rms_emit_y_;
    double sigma_p = rms_dp_p_;
    double beta_x = twiss.bet_x;
    double beta_y = twiss.bet_y;
    double dx = twiss.disp_x;
    double dy = twiss.disp_y;
    double dpx = twiss.disp_dx;
    double dpy = twiss.disp_dy;

    double y_amp = sqrt(2.0*emit_y*beta_y);
    double yp_amp = sqrt(2.0*emit_y/beta_y);
    double x_amp = sqrt(2.0*emit_x*beta_x);
    double xp_amp = sqrt(2.0*emit_x/beta_x);

    double ds_amp, dp_amp;
    if(bunched()){  //bunched beam
        ds_amp = sqrt(2.0)*rms_sigma_s_;
        dp_amp = sqrt(2.0)*sigma_p;
    }

    int cnt = 0;
    for(int i=0; i<n_tr_; ++i){
        double y_spl_tmp = y_amp*y_spl[i];
        double yp_spl_tmp = yp_amp*yp_spl[i];
        for(int j=0; j<n_tr_; ++j){
            double x_spl_tmp = x_amp*x_spl[j];
            double xp_spl_tmp = xp_amp*xp_spl[j];
            if(bunched()) {  //bunched beam
                for(int k=0; k<n_l_; ++k){
                    double ds_spl_tmp = ds_amp*ds_spl[k];
                    double dp_spl_tmp = dp_amp*dp_p_spl[k];
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    ds[cnt] = ds_spl_tmp;
                    dp_p[cnt] = dp_spl_tmp;
                    ++cnt;
                }
            }
            else{   //coasting beam, ds=s-s0 is set to be zero!
                for(int k=-1; k<3; k +=2){
                    x_bet[cnt] = x_spl_tmp;
                    xp_bet[cnt] = xp_spl_tmp;
                    y_bet[cnt] = y_spl_tmp;
                    yp_bet[cnt] = yp_spl_tmp;
                    dp_p[cnt] = k*sigma_p;
                    ++cnt;
                }
            }
        }
    }

    ::adjust_disp(dx, x_bet, dp_p, x, cnt);
    ::adjust_disp(dy, y_bet, dp_p, y, cnt);
    ::adjust_disp(dpx, xp_bet, dp_p, xp, cnt);
    ::adjust_disp(dpy, yp_bet, dp_p, yp, cnt);

    if(bunched())
        update_bet_s();
}
