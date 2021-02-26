#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include <vector>

class Beam;
class IonBeam;
class Ring;
class ElectronBeam;
class Cooler;
class IBSSolver;
class ECoolRate;
class LuminositySolver;

using std::vector;

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
    bool edge_effect = false;
    bool fixed_bunch_length = false;
    virtual void update_ibeam(IonBeam& ionBeam, ElectronBeam& ebeam, double dt)=0;
    virtual void adjust_rf_voltage() = 0;
    virtual void save_ions(double t, IonBeam& ion_sample) { };
    virtual double calc_timestep(double time, int n_steps) const;
    
    DataSink *datasink = nullptr;
    
    Ring &ring;
    Cooler &cooler;

    IBSSolver *ibs_solver = nullptr;
    ECoolRate *ecool_solver = nullptr;
    LuminositySolver *lum_solver = nullptr;
    
    State state;
public:
    Simulator(Ring &_ring, Cooler &_cooler, IBSSolver *_ibs_solver, ECoolRate *_ecool_solver, LuminositySolver *_lum_solver)
    : ring(_ring), cooler(_cooler), ibs_solver(_ibs_solver), ecool_solver(_ecool_solver), lum_solver(_lum_solver) { }
    void set_edge_effect(bool b){edge_effect = b;}
    void set_datasink(DataSink *_datasink) { datasink = _datasink; }

    void set_fixed_bunch_length(bool b){fixed_bunch_length = b; }
    virtual void run(IonBeam& ionBeam,
                     ElectronBeam& ebeam,
                     double time,
                     int n_steps = 0,
                     int state_output_interval = 1,
                     int ion_output_interval = 1
                    );
};

#endif // SIMULATOR_H_INCLUDED
