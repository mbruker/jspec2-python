#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <thread>

class Beam;
class Ions;
class Ring;
class EBeam;
class Cooler;
class IBSSolver;
class ECoolRate;
class FrictionForceSolver;
class LuminositySolver;

using std::string;
using std::ofstream;
using std::vector;

enum class DynamicModel {RMS, PARTICLE, MODEL_BEAM = PARTICLE, TURN_BY_TURN};

//extern double vl_emit_nx, vl_emit_ny, vl_dp_p, vl_sigma_s, vl_rx_ibs, vl_ry_ibs, vl_rs_ibs,
//    vl_rx_ecool, vl_ry_ecool, vl_rs_ecool, vl_rx_total, vl_ry_total, vl_rs_total, vl_t;

class DataSink;

class Simulator{
public:
    struct State {
        double t=0;
        double ex=0, ey=0, dp_p=0, sigma_s=0;
        double rx=0, ry=0, rs=0;
        double rx_ibs=0, ry_ibs=0, rs_ibs=0;
        double rx_ecool=0, ry_ecool=0, rs_ecool=0;
        double rf_voltage=0;
        double luminosity=0;
    };
protected:
//    double t = 0;
    bool edge_effect = false;
    bool fixed_bunch_length = false;
//    bool reset_time = true;
//    bool overwrite = false;
//    bool calc_luminosity = false;
//    int output_itvl = 1;
//    int ion_save_itvl = -1;
//    ofstream outfile;
//    string outfilename = "output_dynamic.txt";
//    vector<double> r_ibs = {0,0,0};
//    vector<double> r_ecool = {0,0,0};
//    vector<double> r = {0,0,0};
//    vector<double> emit = {0,0,0,0};
//    void output_sddshead();
//    void output_to_file();
//    void output(bool bunched=true, double v_rf=0, double lum=0);
    virtual void update_ibeam(Beam& ion, Ions& ion_sample, EBeam& ebeam, double dt)=0;
    virtual void adjust_rf_voltage() = 0;
//    virtual void save_ions(int i, Ions& ion_sample) = 0;
    virtual double calc_timestep(double time, int n_steps) const;
    
    DataSink *datasink = nullptr;
    
    Ring &ring;
    Cooler &cooler;

    IBSSolver *ibs_solver = nullptr;
    ECoolRate *ecool_solver = nullptr;
    FrictionForceSolver *force_solver = nullptr;
    LuminositySolver *lum_solver = nullptr;
    
    State state;
public:
    Simulator(Ring &_ring, Cooler &_cooler, IBSSolver *_ibs_solver, ECoolRate *_ecool_solver, FrictionForceSolver *_force_solver, LuminositySolver *_lum_solver)
    : ring(_ring), cooler(_cooler), ibs_solver(_ibs_solver), ecool_solver(_ecool_solver), force_solver(_force_solver), lum_solver(_lum_solver) { }
    void set_edge_effect(bool b){edge_effect = b;}
//    void set_ion_save(int x){ion_save_itvl = x; }
//    void set_output_file(string filename){outfilename = filename; }
//    void set_output_intvl(int x){output_itvl = x; }

// XXX
void set_datasink(DataSink *_datasink) { datasink = _datasink; }

    void set_fixed_bunch_length(bool b){fixed_bunch_length = b; }
//    void set_ini_time(double t){t0 = t;}
//    void set_reset_time(bool b){reset_time = b;}
//    void set_overwrite(bool b) {overwrite = b; }
//    void set_calc_lum(bool b) {calc_luminosity = b; }

    virtual void precondition(Ions& ion_sample){};
    virtual void run(Beam& ion,
                     Ions& ion_sample,
                     EBeam& ebeam,
                     double time,
                     int n_steps = 0,
                     int state_output_interval = 1,
                     int ion_output_interval = 1
                    );
};

class DataSink {
public:
    virtual void output_simulator_state(const Simulator::State &state) = 0;
};


#endif // SIMULATOR_H_INCLUDED
