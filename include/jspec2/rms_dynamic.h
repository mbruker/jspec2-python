#ifndef RMS_DYNAMIC_H_INCLUDED
#define RMS_DYNAMIC_H_INCLUDED

#include "jspec2/simulator.h"


class RMSModel:public Simulator {
 private:
    void update_ibeam(IonBeam& ionBeam,  ElectronBeam& ebeam, double dt) override;
    void adjust_rf_voltage() override;
 public:
    using Simulator::Simulator;
};



#endif // RMS_DYNAMIC_H_INCLUDED
