#ifndef IONS_HPP
#define IONS_HPP

#include <string>
#include <vector>
#include <array>

#include "jspec2/beam.h"
#include "jspec2/twiss.h"

class Cooler;

using std::vector;
using std::string;
using std::array;

class Ions{
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
public:
    Ions(const Twiss &_twiss)
        : twiss(_twiss) { }
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
    
    const Twiss& get_twiss() const {return twiss;}
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
    virtual void create_samples(const Beam& ion) = 0;
    double beta_s() const {return beta_s_;}
    void update_bet_s(const Beam& ion){beta_s_ = ion.sigma_s()/ion.dp_p();}
};

class Ions_MonteCarlo: public Ions{
private:
    double statistical_emittance_tr(const vector<double>& x, const vector<double>&xp) const;
    double statistical_emittance_l(const vector<double>& dp_p) const;
public:
    Ions_MonteCarlo(const Twiss &_twiss, int n);
    Ions_MonteCarlo(const Twiss &_twiss, std::string filename, int n,int skip = 0, bool binary = false, int n_buffer = 1000);
    virtual array<double,3> emit() const override;
    virtual array<double,3> emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>&ds) const override;
    virtual void create_samples(const Beam& ion);

};

class Ions_SingleParticle: public Ions{
private:
    int n_tr_ = 0;
    int n_l_ = 0;
    vector<double> x_spl, xp_spl, y_spl, yp_spl, ds_spl, dp_p_spl;
public:
    void single_particle_grid(const Beam& ion);
    Ions_SingleParticle(const Twiss &_twiss, int n_tr, int n_l);
    virtual array<double,3> emit() const override;
    virtual array<double,3> emit(const vector<double>& x_bet, const vector<double>& xp_bet, const vector<double>& y_bet, const vector<double>& yp_bet,
                      const vector<double>& dp_p, const vector<double>&ds) const override;
    virtual void create_samples(const Beam& ion);
};

void adjust_disp(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
void adjust_disp_inv(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
#endif
