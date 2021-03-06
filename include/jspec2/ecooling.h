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
class DataSink;

using std::vector;
using std::initializer_list;

class ECoolRate{
protected:
    double bunch_separate_ = 0;
    double t_cooler_ = 0;
    int n_long_sample_ = 50;
    vector<double> ne;
    vector<double> xp_bet, yp_bet, xp, yp, dp_p, x, y, x_bet, y_bet;
    vector<double> v_tr, v_long;
    vector<double> force_x, force_y, force_z;

    FrictionForceSolver *force_solver;
    FrictionForceSolver *longitudinal_force_solver;
    DataSink *datasink = nullptr;

    void electron_density(const IonBeam& ion_sample, ElectronBeam &ebeam);
    void init_scratch(int n_sample);
    void space_to_dynamic(const IonBeam &ion_sample);
    void beam_frame(double gamma_e);
    void lab_frame(double gamma_e);
    void force(const IonBeam &ion, const ElectronBeam &ebeam, const Cooler &cooler);
//    void restore_velocity(ElectronBeam &ebeam);
    void bunched_to_coasting(IonBeam &ion, ElectronBeam &ebeam, const Cooler &cooler);
    void force_distribute(const IonBeam &ion);
    void apply_kick(const IonBeam& ion);
public:
    ECoolRate(FrictionForceSolver *_force_solver, FrictionForceSolver *_longitudinal_force_solver = nullptr)
        : force_solver(_force_solver), longitudinal_force_solver(_longitudinal_force_solver) { }
    void set_force_datasink(DataSink *_datasink) { datasink = _datasink; }
    void adjust_rate(const IonBeam &ion, const ElectronBeam &ebeam, initializer_list<double*> func);
    
    // getters for computation vectors
    const vector<double> &get_force_x() const { return force_x; }
    const vector<double> &get_force_y() const { return force_y; }
    const vector<double> &get_force_z() const { return force_z; }
    
    double t_cooler() const {return t_cooler_;}
    void set_n_long_sample(int n){n_long_sample_ = n;}
    rate3d ecool_rate(IonBeam &ion, const Cooler &cooler, ElectronBeam &ebeam, const Ring &ring);
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
    void output_force(const IonBeam &ion, const Cooler &cooler, ElectronBeam &ebeam);
//    ForceCurve():ECoolRate(){save_force = true;}
};

#endif // ECOOLING_H
