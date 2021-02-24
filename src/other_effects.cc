#include <vector>

#include "jspec2/other_effects.h"
#include "jspec2/electron_beam.h"
#include "jspec2/ion_beam.h"
#include "jspec2/cooler.h"

using std::vector;

// TODO This function does not seem to be needed, and it doesn't do anything.
void edge_effect(ElectronBeam& ebeam, IonBeam& ion, Cooler& cooler, double dt) {
    const int n_sample = ion.n_sample();
    vector<double> ez(n_sample);
    ebeam.edge_field(cooler, ion.cdnt_x(), ion.cdnt_y(), ion.cdnt_ds(), ez, n_sample);
    vector<double>& dp_p = ion.cdnt_dp_p();
    double p0 = ion.p0_SI();
    double inv_p0 = 1/p0;
    for(int i=0; i<n_sample; ++i) {
        dp_p.at(i) = (dp_p.at(i)*p0 + ez.at(i)*dt)*inv_p0;
    }
}
