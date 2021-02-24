#ifndef TURN_BY_TURN_H_INCLUDED
#define TURN_BY_TURN_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>

#include "jspec2/particle_model.h"

class TurnByTurnModel: public ParticleModel {
protected:
    void move_particles(IonBeam& ionBeam) override;
    void apply_edge_kick(ElectronBeam& ebeam, IonBeam& ionBeam) override;
    virtual double calc_timestep(double time, int n_steps) const override;

public:
    using ParticleModel::ParticleModel;
};

#endif // TURN_BY_TURN_H_INCLUDED
