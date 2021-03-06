#ifndef BEAM_H
#define BEAM_H

#include <cstdio>
#include <string>
#include <vector>

#include "jspec2/constants.h"
#include "jspec2/arbitrary_electron_beam.h"

class Cooler;
using std::vector;


class ElectronBeam {
public:
    enum class Shape {UNIFORM_CYLINDER, GAUSSIAN_BUNCH, UNIFORM_BUNCH, GAUSSIAN_CYLINDER, ELLIPTIC_UNIFORM_BUNCH,
        UNIFORM_HOLLOW, UNIFORM_HOLLOW_BUNCH, PARTICLE_BUNCH};
    
    enum class Velocity {CONST, USER_DEFINE, SPACE_CHARGE, VARY, VARY_X, VARY_Y, VARY_Z}  ;
    enum class Temperature {CONST, USER_DEFINE, SPACE_CHARGE, VARY, VARY_X, VARY_Y, VARY_Z}  ;
    enum class EdgeEffect {Rising, Falling};
protected:
    double kinetic_energy_ = 0;
    double gamma_ = 1;
    double beta_ = 0;
    bool bunched_ = true;
    double center_x_ = 0;
    double center_y_ = 0;
    double center_z_ = 0;
    double neutralisation_ = 2;
    Velocity velocity_ = Velocity::CONST;
    Temperature temperature_ = Temperature::CONST;
    vector<double> tpr_t;
    vector<double> tpr_l;
    vector<double> v_rms_t;
    vector<double> v_rms_l;
    vector<double> v_avg_x;
    vector<double> v_avg_y;
    vector<double> v_avg_l;
    bool multi_bunches_ = false;
    int n_; //Number of bunches
    vector<double> cx_;     //List of cxs.
    vector<double> cy_;     //List of cys.
    vector<double> cz_;     //List of czs.
    bool p_shift_ = false;             //Position shift. false: ion center and e- center overlap, true: there's a shift between the beam
    bool v_shift_ = false;             //Velocity shift.
    double cv_l_ = 0;
public:
    virtual ~ElectronBeam(){};
    Velocity velocity() const {return velocity_;}
    Temperature temperature() const {return temperature_;}
    int charge_number() const {return -1;}
    double mass() const {return k_me;}
    double mass_number() const {return k_me/k_u;}
    double mass_SI() const {return k_me*1e6*k_e;}
    double kinetic_energy() const {return kinetic_energy_;}
    double gamma() const {return gamma_;}
    double beta() const {return beta_;}
    bool bunched()const {return bunched_;}
    void set_p_shift(bool b){p_shift_ = b;}
    void set_v_shift(bool b){v_shift_ = b;}
    bool p_shift() const {return p_shift_;}
    bool v_shift() const {return v_shift_;}
    void set_cv_l(double x){cv_l_ = x; v_shift_ = true;}
    double cv_l() const {return cv_l_;}
    virtual Shape shape() const = 0;
    virtual double length() const = 0;
    double neutral() const {return neutralisation_;}
    
    // convenience functions for the UI. Make sure one of them is called
    void set_kinetic_energy(double ke);
    void set_gamma(double g);
    void set_beta(double b);
    
    void set_center(double cx, double cy, double cz);
    void set_tpr(double tpr_tr, double trp_long);
    void set_v_rms(double v_rms_tr, double v_rms_long);
    void set_v_avg(double v_avg_tx, double v_avg_ty, double v_avg_long);
    void set_neutral(double x){neutralisation_ = x;}
    void set_multi_bunches(bool b){multi_bunches_ = b;}
    bool multi_bunches() const {return multi_bunches_;}
    const vector<double>& get_v_rms_tr() const { return v_rms_t; };
    const vector<double>& get_v_rms_l() const { return v_rms_l; }
    virtual void edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                             vector<double>& field, int n){};
    virtual void edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                             vector<double>& field, int n, double cx, double cy, double cz){};
    void set_n_bunches(int n){n_ = n; cx_.resize(n); cy_.resize(n), cz_.resize(n);}
    virtual void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) = 0;
    virtual void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                         double cx, double cy, double cz) = 0;
    void multi_density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n);
    void multi_density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                         double cx, double cy, double cz);
    void multi_edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                             vector<double>& field, int n);
    void multi_edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                             vector<double>& field, int n, double cx, double cy, double cz);
};

class UniformCylinder: public ElectronBeam{
    double current_;                   //Current of the beam in A
    double radius_;              //Radius of the beam in meter
 public:
    void density (const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density (const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    double current() const {return current_;}
    double radius() const {return radius_;}
    Shape shape() const {return Shape::UNIFORM_CYLINDER;}
    double length() const {perror("length() not defined for UniformCylinder, which is coasting"); return 0;}
    UniformCylinder(double current, double radius, double neutralisation=2):current_(current),radius_(radius)
                    {bunched_ = false;};
};

class UniformHollow: public ElectronBeam {
    double current_;    //Peak current, the current as if the beam is coasting.
    double in_radius_;
    double out_radius_;
 public:
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    double current() const {return current_;}
    double out_radius() const {return out_radius_;}
    double in_radius() const {return in_radius_;}
    Shape shape() const {return Shape::UNIFORM_HOLLOW;}
    double length() const {perror("length() not defined for UniformHollow, which is coasting"); return 0;}
    UniformHollow(double current, double in_radius, double out_radius):current_(current),
        in_radius_(in_radius), out_radius_(out_radius){bunched_ = false;};
};


class UniformHollowBunch: public ElectronBeam {
    double current_;
    double in_radius_;
    double out_radius_;
    double length_;
 public:
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    double current() const {return current_;}
    double out_radius() const {return out_radius_;}
    double in_radius() const {return in_radius_;}
    Shape shape() const {return Shape::UNIFORM_HOLLOW_BUNCH;}
    double length() const {return length_;}
    UniformHollowBunch(double current, double in_radius, double out_radius, double length):current_(current),
        in_radius_(in_radius), out_radius_(out_radius), length_(length) {}
};


class UniformBunch: public ElectronBeam{
    double current_;                   //Current of the beam in A, assuming the beam is DC.
    double radius_;              //Radius of the beam in meter
    double length_;
    double t_rising_ = 0;
    double t_falling_ = 0;
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    Shape shape() const {return Shape::UNIFORM_BUNCH;}
    double length() const {return length_;}
    double current() const {return current_;}
    double radius() const {return radius_;}
    void set_rising_time(double x){t_rising_ = x;}
    void set_falling_time(double x){t_falling_ = x;}
    void edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                    vector<double>& field, int n);
    void edge_field(const Cooler& cooler, const vector<double>&x, const vector<double>& y, const vector<double>&z,
                    vector<double>& field, int n, double cx, double cy, double cz);
    UniformBunch(double current, double radius, double length):current_(current),radius_(radius),
            length_(length){};

};


class EllipticUniformBunch: public ElectronBeam{
    double current_;
    double rh_;         //half horizontal axis
    double rv_;         //half vertical axis
    double length_;     //bunch length
public:
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    Shape shape() const {return Shape::ELLIPTIC_UNIFORM_BUNCH;}
    double length() const {return length_;}
    EllipticUniformBunch(double current, double rh, double rv, double length):current_(current),
            rh_(rh),rv_(rv),length_(length){};
};


class GaussianBunch: public ElectronBeam{
    double n_electron_;
    double sigma_x_;
    double sigma_y_;
    double sigma_s_;
    double sigma_xp_;
    double sigma_yp_;
    double sigma_dpp_;
 public:
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    Shape shape() const {return Shape::GAUSSIAN_BUNCH;}
    double length() const {return 6*sigma_s_;}
    void set_angles(double sigma_xp, double sigma_yp, double sigma_dpp);
    GaussianBunch(double n_electron, double sigma_x, double sigma_y, double sigma_s):n_electron_(n_electron),
                sigma_x_(sigma_x),sigma_y_(sigma_y),sigma_s_(sigma_s){};


};

class ParticleBunch: public ElectronBeam {
    double n_electron_;
    std::string filename_;
    long int n_ = 0;
    double length_ = 0;
    bool v_x_corr_ = false;    //Velocity position correlation
    int line_skip_ = 0;
    vector<Box> tree_;
    vector<long int> list_e_;
    int s_ = 100;
    bool binary_ = false;
    int buffer_ = 1000;
public:
    std::vector<double> x, y, z, vx, vy, vz;  //Electron phase space coordinates
    //Calculate the charge density for a given position (x,y,z) in Lab frame.
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n) override;
    void density(const vector<double>& x, const vector<double>& y, const vector<double>& z, vector<double>& ne, int n,
                double cx, double cy, double cz) override;
    Shape shape() const {return Shape::PARTICLE_BUNCH;}
    double length() const {return length_;}
    bool bunched() const {return true;}
    bool corr() const {return v_x_corr_;}
    void set_corr(bool corr = true){v_x_corr_ = corr;}
    void set_buffer(int n) {buffer_ = n;}
    void set_s(int s) {s_ = s;}
    void set_binary(bool b) {binary_ = b;}
    void set_skip(int n) {line_skip_ = n;}

    // TODO Why doesn't this constructor load the data?
    ParticleBunch(double n_electron, std::string filename, double length = 0)
        : n_electron_(n_electron),
          filename_(filename),
          length_(length)
    {
        temperature_ = Temperature::VARY;
    }
    void load_particle(long int n);
    void load_particle();

};

//class MultiBunches: public ElectronBeam{
//    int n_; //Number of bunches
//    vector<double> cx_;     //List of cxs.
//    vector<double> cy_;     //List of cys.
//    vector<double> cz_;     //List of czs.
//
// public:
//    ElectronBeam* bunches_;
////    ElectronBeam* bunch(){return bunches_;}
//    vector<double>& cx(){return cx_;}
//    vector<double>& cy(){return cy_;}
//    vector<double>& cz(){return cz_;}
//    MultiBunches(int n):n_(n){cx_.resize(n); cy_.resize(n); cz_.resize(n);}
//    ~MultiBunches(){delete bunches_;}
//    MultiBunches(const MultiBunches& obj) = delete;
//    MultiBunches& operator=(const MultiBunches& obj) = delete;
//    Shape shape() const {return bunches_->shape();}
//    double length() const {return bunches_->length();}
//    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n);
//    void density (vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& ne, int n,
//                double cx, double cy, double cz);
//};

#endif // BEAM_H
