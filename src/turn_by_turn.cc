#include <chrono>
#include <cmath>
#include "jspec2/turn_by_turn.h"
#include "jspec2/electron_beam.h"
#include "jspec2/constants.h"
#include "jspec2/ecooling.h"
#include "jspec2/functions.h"
#include "jspec2/other_effects.h"
#include "jspec2/particle_model.h"
#include "jspec2/ring.h"
#include "jspec2/ion_beam.h"
#include "jspec2/cooler.h"

void TurnByTurnModel::apply_edge_kick(ElectronBeam& ebeam, IonBeam& ionBeam) {
//    ::edge_effect(ebeam, ionBeam, cooler, dt);
    const vector<double>& x = ionBeam.cdnt_x();
    const vector<double>& y = ionBeam.cdnt_y();
    const vector<double>& ds = ionBeam.cdnt_ds();
    vector<double>& dp_p = ionBeam.cdnt_dp_p();
    const double p0 = ionBeam.p0_SI();
    const int n = ionBeam.n_sample();
    vector<double> field(n);

//    if(ecool_solver.p_shift_) {
//        if(ebeam.multi_bunches())
//            ebeam.multi_edge_field(x, y, ds, field, n, ion_sample.center_x(), ion_sample.center_y(), ion_sample.center_z());
//        else
//            ebeam.edge_field(x, y, ds, field, n, ion_sample.center_x(), ion_sample.center_y(), ion_sample.center_z());
//    }
//    else {
//        if(ebeam.multi_bunches()) ebeam.multi_edge_field(x, y, ds, field, n);
//        else ebeam.edge_field(cooler, x, y, z, field, n);
//    }

    ebeam.edge_field(cooler, x, y, ds, field, n);
    const double q = ionBeam.charge_number()*k_e;
    const double coef = q*ecool_solver->t_cooler()*cooler.section_number()/p0;
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for(int i=0; i<n; ++i) {
        dp_p[i] *= (1+field.at(i)*coef/dp_p.at(i));
//        dp_p.at(i) = dp_p.at(i)*exp(field.at(i)*coef/dp_p.at(i));

//        double dp = field.at(i)*coef/dp_p.at(i);
//        dp_p[i] = dp>0.15?dp_p[i]*(1+dp):dp_p[i]*exp(dp);
    }
}

void TurnByTurnModel::move_particles(IonBeam& ionBeam) {
    //Transverse
    //New betatron oscillation coordinates
    const Twiss& twiss = ionBeam.get_twiss();

    //Transverse motion by tunes
    const double Qx = ring.qx();
    const double Qy = ring.qy();
    assert(Qx>0 && Qy>0 &&"Transverse tunes are needed for Turn_by_turn model");
    
    vector<double>& x_bet = ionBeam.cdnt_x_bet();
    vector<double>& xp_bet = ionBeam.cdnt_xp_bet();
    vector<double>& y_bet = ionBeam.cdnt_y_bet();
    vector<double>& yp_bet = ionBeam.cdnt_yp_bet();

    vector<double>& dp_p = ionBeam.cdnt_dp_p();
    vector<double>& ds = ionBeam.cdnt_ds();

    int n_sample = ionBeam.n_sample();
    ionBeam.adjust_disp_inv();

    double bet_x = twiss.bet_x;
    double bet_y = twiss.bet_y;
    double alf_x = twiss.alf_x;
    double alf_y = twiss.alf_y;
    double gamma_x = (1+alf_x*alf_x)/bet_x;
    double gamma_y = (1+alf_y*alf_y)/bet_y;
    const double sin_phi_x = sin(2*k_pi*Qx);
    const double cos_phi_x = cos(2*k_pi*Qx);
    const double sin_phi_y = sin(2*k_pi*Qy);
    const double cos_phi_y = cos(2*k_pi*Qy);
    
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif // _OPENMP
    for (int i=0; i<n_sample; ++i) {
        const double x_bet_0 = x_bet[i];
        const double xp_bet_0 = xp_bet[i];
        x_bet[i] = (cos_phi_x+alf_x*sin_phi_x)*x_bet_0 + bet_x*sin_phi_x*xp_bet_0;
        xp_bet[i] = -gamma_x*sin_phi_x*x_bet_0 + (cos_phi_x-alf_x*sin_phi_x)*xp_bet_0;
        const double y_bet_0 = y_bet[i];
        const double yp_bet_0 = yp_bet[i];
        y_bet[i] = (cos_phi_y+alf_y*sin_phi_y)*y_bet_0 + bet_y*sin_phi_y*yp_bet_0;
        yp_bet[i] = -gamma_y*sin_phi_y*y_bet_0 + (cos_phi_y-alf_y*sin_phi_y)*yp_bet_0;
    }

    //Longitudinal motion.
//    vector<double>& dp_p = ionBeam.cdnt_dp_p();
//    vector<double>& ds = ionBeam.cdnt_ds();
    if (ring.qs()>0||ring.rf_voltage()>0) {    //RF, synchrotron oscillation.
//        assert(ring.tunes->qs>0||ring.rf->v>0&&"Longitudinal tune or RF cavity needed for Turn_by_turn model");

        if(ring.rf_voltage()>0) { //Longitudinal motion by RF.
            double circ = ring.circ();
            double beta2 = ionBeam.beta()*ionBeam.beta();
            double beta2_inv = 1/beta2;
            double total_energy = ionBeam.gamma()*ionBeam.mass(); //ion total energy [MeV/c^2]
            double total_energy_inv = 1/total_energy;
            double adj_dp2dE = beta2*total_energy;

            double volt = ring.rf_voltage();
            double phi_s = ring.rf_phi();
//            double phi_0 = ring.rf->phi_0();
            double h = ring.rf_h();
//            double s_s = phi_s*circ/(h*2*k_pi);
            double half_phase = h*k_pi;
            double total_phase = h*2*k_pi;
            double adj_s2phi = total_phase/circ;
            double adj_phi2s = 1/adj_s2phi;
            double adj_dE = ionBeam.charge_number()*volt*1e-6; // [MeV/c^2]
    //                double adj_dE = ion.charge_number()*ring.rf_->volt()*1e-6; // [MeV/c^2]
            double sin_phi_s = sin(phi_s);
//            double sin_phi_s = sin(phi_s+phi_0);
            double eta = ring.slip_factor();
            double adj_dE2dphi = total_phase*eta*beta2_inv*total_energy_inv;
            #ifdef _OPENMP
                #pragma omp parallel for
            #endif // _OPENMP
            for(int i = 0; i < n_sample; ++i) {
                dp_p[i] *= adj_dp2dE; //dp/p -> dE/E -> dE in [MeV/c^2]
//                ds[i] += s_s;  //s = ds + s_s: adjust ds to be measured from the start of the ring
                ds[i] *= adj_s2phi;  //phi = s*h*2*pi/circ: s -> phi
//                dp_p[i] += adj_dE*(sin(ds[i]+phi_0)-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                dp_p[i] += adj_dE*(sin(ds[i])-sin_phi_s); //dE_n+1 = dE_n + q*v*[sin(phi)-sin(phi_s)]
                ds[i] += adj_dE2dphi*dp_p[i]; //phi_n+1 = phi_n + 2*pi*h*eta*dE_n+1 / (beta^2 * E)
                ds[i] += half_phase;    //phi in [0, total phase]
                ds[i] = fmod(ds[i], total_phase);  //phi_n+1 module h*2*pi
                if (ds[i]<0) ds[i] += total_phase; //adjust phi is phi is less than zero
                ds[i] -= half_phase;    //phi in [-half_phase, half_phase]
                ds[i] *= adj_phi2s;  //phi -> s
//                ds[i] -= s_s;        // ds = s - s_s: adjust ds back to be centered about s_s
                dp_p[i]  = dp_p[i]*total_energy_inv*beta2_inv; //dE -> dE/E -> dp/p = beta*beta*dE/E;
            }
        }
        else if(ring.qs()>0) {//Longitudinal motion by tune
            const double sin_phi = sin(2*k_pi*ring.qs());
            const double cos_phi = cos(2*k_pi*ring.qs());
            const double beta_s = ionBeam.beta_s();
            const double inv_beta_s = 1/beta_s;
            #ifdef _OPENMP
                #pragma omp parallel for
            #endif // _OPENMP
            for (int i=0; i<n_sample; ++i) {
                const double dp_p_0 = dp_p[i];
                const double ds_0 = ds[i];
                dp_p[i] = cos_phi*dp_p_0 - sin_phi*ds_0*inv_beta_s;
                ds[i] = sin_phi*dp_p_0*beta_s + cos_phi*ds_0;
            }
        }
    }
    else {  //No RF.
        double gamma_0 = ionBeam.gamma();
        double beta_0 = ionBeam.beta();
        double half_length = 0.5*ring.circ();
        #ifdef _OPENMP
            #pragma omp parallel for
        #endif // _OPENMP
        for(int i=0; i<n_sample; ++i) {
            double gamma2 = 1+(1+dp_p[i])*(1+dp_p[i])*(gamma_0*gamma_0-1);
            double beta = sqrt(1-1/gamma2);
            double s = (beta/beta_0-1)*2*half_length;
            ds[i] += half_length;    //s in [0, 2*half_length]
            ds[i] += s;
            ds[i] = fmod(ds[i], 2*half_length);
            if(ds[i]<0) ds[i] += 2*half_length;
            ds[i] -= half_length;   //s in [-half_length, half_length]
        }
    }

    ionBeam.adjust_disp();
}

double TurnByTurnModel::calc_timestep(double /*time*/, int /*n_steps*/) const
{
    return 1.0 / ring.f0();
}
