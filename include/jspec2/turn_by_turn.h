#ifndef TURN_BY_TURN_H_INCLUDED
#define TURN_BY_TURN_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>

#include "jspec2/particle_model.h"

class TurnByTurnModel: public ParticleModel {
protected:
    void move_particles(Beam& ion, Ions& ion_sample) override;
    void apply_edge_kick(EBeam& ebeam, Beam& ion, Ions& ion_sample) override;
    ofstream out_single_particle;
    string filename_single_particle = "single_particle_";
    int idx = -1;
    virtual double calc_timestep(double time, int n_steps) const override;

public:
    using ParticleModel::ParticleModel;
};

#endif // TURN_BY_TURN_H_INCLUDED
