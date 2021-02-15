#ifndef RMS_DYNAMIC_H_INCLUDED
#define RMS_DYNAMIC_H_INCLUDED

#include "jspec2/simulator.h"


class RMSModel:public Simulator {
 private:
    void update_ibeam(Beam& ion, Ions& ion_sample,  EBeam& ebeam, double dt) override;
    void adjust_rf_voltage() override;
//    void save_ions(int i, Ions& ion_sample) override { };
 public:
    using Simulator::Simulator;
};



#endif // RMS_DYNAMIC_H_INCLUDED
