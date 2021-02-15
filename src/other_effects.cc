#include <vector>

#include "jspec2/other_effects.h"
#include "jspec2/beam.h"
#include "jspec2/ions.h"
#include "jspec2/cooler.h"

using std::vector;

void edge_effect(EBeam& ebeam, Beam& ion, Ions& ion_sample, Cooler& cooler, double dt) {
    int n_sample = ion_sample.n_sample();
    vector<double> ez(n_sample);
    ebeam.edge_field(cooler, ion_sample.cdnt_x(), ion_sample.cdnt_y(), ion_sample.cdnt_ds(), ez, n_sample);
    vector<double>& dp_p = ion_sample.cdnt_dp_p();
    double p0 = ion.p0_SI();
    double inv_p0 = 1/p0;
    for(int i=0; i<n_sample; ++i) {
        dp_p.at(i) = (dp_p.at(i)*p0 + ez.at(i)*dt)*inv_p0;
    }
}
