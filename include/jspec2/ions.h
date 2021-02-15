#ifndef IONS_HPP
#define IONS_HPP

#include <string>
#include <vector>

#include "jspec2/ring.h"
class Cooler;

using std::vector;
using std::string;

enum class IonSampleType {MONTE_CARLO, USER_DEFINE, SINGLE_PARTICLE};
enum class Phase {X_BET, XP_BET, Y_BET, YP_BET, X, Y, XP, YP, DS, DP_P};

class Ions{
protected:
    vector<double> x_bet, xp_bet, y_bet, yp_bet;
    vector<double> x, y, xp, yp, ds, dp_p;
    int n_ = 0; //Number of sample particles.
    IonSampleType sample_type_ = IonSampleType::MONTE_CARLO;
    Twiss twiss;
    double center_[3] = {0,0,0};
    bool bunched_ = true;
    double beta_s_ = 0;
public:
    void adjust_disp();
    void adjust_disp_inv();
    void save_ions_sdds(string filename) const;
    const vector<double>& get_cdnt(Phase p) const;
    vector<double>& cdnt(Phase p);
    const Twiss& get_twiss() const {return twiss;}
//    void set_n_sample(int n){n_ = n;}
    int n_sample() const { return n_; }
    IonSampleType sample_type() const {return sample_type_;}
    void set_twiss(const Twiss& t);
    void set_twiss(const Cooler& cooler);
    void center(double &cx, double &cy, double &cz){cx = center_[0]; cy = center_[1]; cz = center_[2];}
    virtual void emit(double& emit_x, double& emit_y, double& emit_s) = 0;
    virtual void emit(vector<double>& x_bet, vector<double>& xp_bet, vector<double>& y_bet, vector<double>& yp_bet,
                      vector<double>& dp_p,vector<double>&ds, double& emit_x, double& emit_y, double& emit_s) = 0;
    virtual void create_samples(const Beam& ion) = 0;
    double beta_s() const {return beta_s_;}
    void update_bet_s(const Beam& ion){beta_s_ = ion.sigma_s()/ion.dp_p();}
};

class Ions_MonteCarlo: public Ions{
private:
//    double emit(vector<double>& x, vector<double>&xp, int n);
//    double emit_p(vector<double>& dp_p, int n);
public:
    Ions_MonteCarlo(int n);
    Ions_MonteCarlo(std::string filename, int n,int skip = 0, bool binary = false, int n_buffer = 1000);
    virtual void emit(double& emit_x, double& emit_y, double& emit_s);
    virtual void emit(vector<double>& x_bet, vector<double>& xp_bet, vector<double>& y_bet, vector<double>& yp_bet,
                      vector<double>& dp_p, vector<double>&ds, double& emit_x, double& emit_y, double& emit_s);
    virtual void create_samples(const Beam& ion);

};

class Ions_SingleParticle: public Ions{
private:
    int n_tr_ = 0;
    int n_l_ = 0;
    vector<double> x_spl, xp_spl, y_spl, yp_spl, ds_spl, dp_p_spl;
public:
    void single_particle_grid(const Beam& ion);
    Ions_SingleParticle(int n_tr, int n_l);
    virtual void emit(double& emit_x, double& emit_y, double& emit_s);
    virtual void emit(vector<double>& x_bet, vector<double>& xp_bet, vector<double>& y_bet, vector<double>& yp_bet,
                      vector<double>& dp_p, vector<double>&ds, double& emit_x, double& emit_y, double& emit_s);
    virtual void create_samples(const Beam& ion);
};

double emit_p(vector<double>& dp_p, int n);
double emit(vector<double>& x, vector<double>&xp, int n);
void adjust_disp(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
void adjust_disp_inv(double dx, vector<double>& x_bet, vector<double>& dp_p, vector<double>& x, int n);
#endif