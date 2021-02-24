#ifndef OTHER_EFFECTS_H
#define OTHER_EFFECTS_H

// XXX This file smells dirty. Check later

class ElectronBeam;
class Beam;
class Ions;
class Cooler;

void edge_effect(ElectronBeam& ebeam,Beam& ion, Ions& ion_sample, Cooler& cooler, double dt);
#endif // OTHER_EFFECTS_H
