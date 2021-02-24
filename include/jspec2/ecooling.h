#ifndef ECOOLING_H
#define ECOOLING_H

#include <fstream>
#include <initializer_list>
#include <string>
#include <vector>

#include "jspec2/rate.h"

class ElectronBeam;
class IonBeam;
class FrictionForceSolver;
class Cooler;
class Ring;

using std::vector;
using std::initializer_list;

class ECoolRate{
protected:
    double bunch_separate_ = 0;
    double t_cooler_ = 0;
    int n_long_sample_ = 50;
    int scratch_size = 0;
    bool dual_force_solver = false;
    bool save_force = false;
    vector<double> ne;
    vector<double> xp_bet, yp_bet, xp, yp, dp_p, x, y, x_bet, y_bet;
    vector<double> v_tr, v_long;
    vector<double> force_x, force_y, force_z;
    void electron_density(const IonBeam& ion_sample, ElectronBeam &ebeam);
    void init_scratch(int n_sample);
    void space_to_dynamic(int n_sample, const IonBeam &ion_sample);
    void beam_frame(int n_sample, double gamma_e);
    void force(int n_sample, const IonBeam &ion, const ElectronBeam &ebeam, const Cooler &cooler, FrictionForceSolver &force_solver);
    void restore_velocity(int n_sample, ElectronBeam &ebeam);
    void bunched_to_coasting(IonBeam &ion, ElectronBeam &ebeam, const Cooler &cooler, FrictionForceSolver &force_solver);
    void lab_frame(int n_sample, double gamma_e);
    void force_distribute(int n_sample, const IonBeam &ion);
    void apply_kick(int n_sample, const IonBeam& ion);
    void save_force_sdds_head(std::ofstream& of, int n);
    FrictionForceSolver* force_solver_l;
public:
    void set_save_force(bool b){save_force = b;}
    void set_dual_force_solver(bool b){dual_force_solver = b;}
    void set_second_force_solver(FrictionForceSolver* force) {force_solver_l = force;}
    void adjust_rate(const IonBeam &ion, const ElectronBeam &ebeam, initializer_list<double*> func);
    
    // getters for computation vectors
    const vector<double> &get_force_x() const { return force_x; }
    const vector<double> &get_force_y() const { return force_y; }
    const vector<double> &get_force_z() const { return force_z; }
    
    double t_cooler() const {return t_cooler_;}
    void set_n_long_sample(int n){n_long_sample_ = n;}
    rate3d ecool_rate(FrictionForceSolver &force, IonBeam &ion, const Cooler &cooler, ElectronBeam &ebeam, const Ring &ring);
};

class ForceCurve: public ECoolRate {
    int n_tr = 1;
    int n_l = 1;
    double dp_p = 0;    //Longitudinal momentum / reference momentum
    double angle = 0;   //Angle in [rad]
    double density_e = 0;
public:
    void set_n_tr(int n){n_tr = n;}
    void set_n_l(int n){n_l = n;}
    void set_electron_density(double x){density_e = x;}
    void set_dp_p(double x) {dp_p = x;}
    void set_angle(double x) {angle = x;}
    void force_to_file(FrictionForceSolver &force_solver, const IonBeam &ion, const Cooler &cooler, ElectronBeam &ebeam);
//    ForceCurve():ECoolRate(){save_force = true;}
};

#endif // ECOOLING_H
