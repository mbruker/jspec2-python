#ifndef IONS_HPP
#define IONS_HPP

#include <string>
#include <vector>
#include <array>

#include "jspec2/twiss.h"

class Cooler;

using std::vector;
using std::string;
using std::array;

class IonBeam {
protected:
    vector<double> x_bet, xp_bet, y_bet, yp_bet;
    vector<double> x, y, xp, yp, ds, dp_p;
    int n_ = 0; //Number of sample particles.
    Twiss twiss;
    double center_x_ = 0;
    double center_y_ = 0;
    double center_z_ = 0;
    bool bunched_ = true;
    double beta_s_ = 0;

    // kinematic quantities; intended to be immutable
    int charge_number_;   //Number of charges
    double mass_;    //unit in MeV/c^2
    double r_;       //classical radius, in m
    double kinetic_energy_;      //kinetic energy, in MeV
    double beta_;    //Lorentz factors
    double gamma_;   //Lorentz factors
    double particle_number_; //number of particles
    double p0_SI_; //momentum in kg*m/s
    
    // statistical phase space quantities
    double rms_emit_nx_; //normalized horizontal emittance, in m
    double rms_emit_ny_; //normalized vertical emittance, in m
    double rms_emit_x_;  //geometrical horizontal emittance, in m
    double rms_emit_y_;  //geometrical vertical emittance, in m
    double rms_dp_p_;     //momentum spread dp/p
    double rms_energy_spread_;       // dE/E
    double rms_dv_v_;     // dv/v
    double rms_sigma_s_; //RMS bunch length. set it to -1 for coasting beam, in m

public:
    IonBeam(const Twiss &_twiss, int charge_number, double mass, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle);
    void adjust_disp();
    void adjust_disp_inv();
    void save_ions_sdds(string filename) const;
    
    // getters for computation vectors
    vector<double>& cdnt_x_bet() { return x_bet; }
    vector<double>& cdnt_xp_bet() { return xp_bet; }
    vector<double>& cdnt_y_bet() { return y_bet; }
    vector<double>& cdnt_yp_bet() { return yp_bet; }
    vector<double>& cdnt_x() { return x; }
    vector<double>& cdnt_y() { return y; }
    vector<double>& cdnt_xp() { return xp; }
    vector<double>& cdnt_yp() { return yp; }
    vector<double>& cdnt_ds() { return ds; }
    vector<double>& cdnt_dp_p() { return dp_p; }
    const vector<double>& cdnt_x_bet() const { return x_bet; }
    const vector<double>& cdnt_xp_bet() const { return xp_bet; }
    const vector<double>& cdnt_y_bet() const { return y_bet; }
    const vector<double>& cdnt_yp_bet() const { return yp_bet; }
    const vector<double>& cdnt_x() const { return x; }
    const vector<double>& cdnt_y() const { return y; }
    const vector<double>& cdnt_xp() const { return xp; }
    const vector<double>& cdnt_yp() const { return yp; }
    const vector<double>& cdnt_ds() const { return ds; }
    const vector<double>& cdnt_dp_p() const { return dp_p; }
    
    // getters for statistical phase space quantities
    double rms_emit_nx() const { return rms_emit_nx_; }
    double rms_emit_ny() const { return rms_emit_ny_; }
    double rms_emit_x() const { return rms_emit_x_; }
    double rms_emit_y() const { return rms_emit_y_; }
    double rms_dp_p() const { return rms_dp_p_; }
//    double energy_spread() const {return rms_energy_spread_;}
//    double velocity_spread() const {return rms_dv_v_;}
    double rms_sigma_s() const { return rms_sigma_s_; }

        // TODO These setters should likely not exist.
        // They are currently needed to update the rms quantities from outside the class.
    void set_emit_nx(double x)
    {
        rms_emit_nx_ = x;
        rms_emit_x_ = rms_emit_nx_/(beta_*gamma_);
    }
    void set_emit_ny(double x)
    {
        rms_emit_ny_ = x;
        rms_emit_y_ = rms_emit_ny_/(beta_*gamma_);
    }
    void set_emit_x(double x)
    {
        rms_emit_x_ = x;
        rms_emit_nx_ = beta_*gamma_*rms_emit_x_;
    }
    void set_emit_y(double x)
    {
        rms_emit_y_ = x;
        rms_emit_ny_ = beta_*gamma_*rms_emit_y_;
    }
    void set_dp_p(double x)
    {
        rms_dp_p_ = x;
        rms_energy_spread_ = beta_*beta_*rms_dp_p_;
        rms_dv_v_ = rms_dp_p_/(gamma_*gamma_);
    }
    void set_sigma_s(double x)
    {
        rms_sigma_s_ = x;
    }
    
    
    // getters for kinematic quantities and misc. properties
    bool bunched() const { return bunched_; }
    int charge_number() const { return charge_number_; }
    double mass() const { return mass_; }
    double kinetic_energy() const { return kinetic_energy_; }
    double beta() const { return beta_; }
    double gamma() const { return gamma_; }
    double p0_SI() const { return p0_SI_; }
    double p0() const { return beta_*gamma_*mass_; }  //Momentum in MeV/c
    double classical_radius() const { return r_; }
    double particle_number() const { return particle_number_; }
    
    const Twiss& get_twiss() const { return twiss; }
    int n_sample() const { return n_; }
    double center_x() const { return center_x_; }
    void set_center_x(double center_x) { center_x_ = center_x; }
    double center_y() const { return center_y_; }
    void set_center_y(double center_y) { center_y_ = center_y; }
    double center_z() const { return center_z_; }
    void set_center_z(double center_z) { center_z_ = center_z; }
    
    virtual array<double,3> emit() const = 0;
    virtual array<double,3> emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>&ds) const = 0;
    virtual void create_samples() = 0;
    double beta_s() const {return beta_s_;}
    // TODO Check if needed
    void update_bet_s(){beta_s_ = rms_sigma_s_ / rms_dp_p_;}
};

class IonBeam_MonteCarlo: public IonBeam {
private:
    double statistical_emittance_tr(const vector<double>& x, const vector<double>&xp) const;
    double statistical_emittance_l(const vector<double>& dp_p) const;
public:
    IonBeam_MonteCarlo(const Twiss &_twiss, int charge_number, double mass, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle, int n);
    IonBeam_MonteCarlo(const Twiss &_twiss, int charge_number, double mass, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle, std::string filename, int n,int skip = 0, bool binary = false, int n_buffer = 1000);
    virtual array<double,3> emit() const override;
    virtual array<double,3> emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>&ds) const override;
    virtual void create_samples() override;
};

class IonBeam_SingleParticle: public IonBeam {
private:
    int n_tr_ = 0;
    int n_l_ = 0;
    vector<double> x_spl, xp_spl, y_spl, yp_spl, ds_spl, dp_p_spl;
public:
    void single_particle_grid();
    IonBeam_SingleParticle(const Twiss &_twiss, int charge_number, double mass, double kinetic_energy, double emit_nx, double emit_ny, double dp_p,
           double sigma_s, double n_particle, int n_tr, int n_l);
    virtual array<double,3> emit() const override;
    virtual array<double,3> emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>&ds) const override;
    virtual void create_samples() override;
};

void adjust_disp(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
void adjust_disp_inv(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
#endif
