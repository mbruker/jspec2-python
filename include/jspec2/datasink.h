#ifndef DATASINK_H
#define DATASINK_H

#include <vector>
#include "simulator.h"

using std::vector;

class IonBeam;

class DataSink {
public:
    virtual void output_simulator_state(const Simulator::State &state) = 0;
    virtual void output_ion_phasespace(double t, const IonBeam &ionBeam) = 0;
    virtual void output_ecool_force(const vector<double> &n_e, const vector<double> &v_tr, const vector<double> &v_l, const vector<double> &force_x, const vector<double> &force_z) = 0;
};

#endif 
