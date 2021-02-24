#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "jspec2/simulator.h"
#include <vector>

class IonBeam;
class Ring;
class ElectronBeam;
class Cooler;

class ParticleModel: public Simulator {
 protected:
    virtual void update_ibeam(IonBeam& ionBeam, ElectronBeam& ebeam, double dt) override;
    void apply_cooling_kick(double freq, IonBeam& ionBeam, double dt);
    void apply_ibs_kick(IonBeam& ionBeam, double dt);
    void ibs_kick(int n_sample, double rate, double twiss, double emit, vector<double>& v, double dt);
    virtual void move_particles(IonBeam& ionBeam);
    virtual void apply_edge_kick(ElectronBeam& ebeam, IonBeam& ionBeam) { };
    void update_beam_parameters(IonBeam& ionBeam);
    virtual void adjust_rf_voltage() override { };
//    virtual void save_ions(int i, IonBeam& ion_sample) override;
 public:
    using Simulator::Simulator;
};

#endif
