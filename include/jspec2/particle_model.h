#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "jspec2/simulator.h"
#include <vector>

class Beam;
class Ions;
class Ring;
class EBeam;
class Cooler;

class ParticleModel: public Simulator {
 protected:
    virtual void update_ibeam(Beam& ion, Ions& ion_sample, EBeam& ebeam, double dt) override;
    void apply_cooling_kick(double freq, Beam& ion, Ions& ion_model, double dt);
    void apply_ibs_kick(Beam& ion, Ions& ion_sample, double dt);
    void ibs_kick(int n_sample, double rate, double twiss, double emit, vector<double>& v, double dt);
    virtual void move_particles(Beam& ion, Ions& ion_sample);
    virtual void apply_edge_kick(EBeam& ebeam, Beam& ion, Ions& ion_sample) { };
    void update_beam_parameters(Beam &ion, Ions& ion_sample);
    virtual void adjust_rf_voltage() override { };
//    virtual void save_ions(int i, Ions& ion_sample) override;
 public:
    using Simulator::Simulator;
};

#endif // IBS_HPP
