#ifndef COOLER_H
#define COOLER_H

#include "twiss.h"

class Cooler{
    double length_;      // in meter
    double section_number_;
    double magnetic_field_;           // in Tesla
    Twiss twiss_;
    double pipe_radius_;
 public:
    double length() const {return length_;}
    double section_number() const {return section_number_;}
    double magnetic_field() const {return magnetic_field_;}
    const Twiss &twiss() const { return twiss_; }
    double pipe_radius() const {return pipe_radius_;}
    Cooler(double length, double section_number, double magnetic_field, const Twiss &twiss, double pipe_radius = 0)
        : length_(length),
          section_number_(section_number),
          magnetic_field_(magnetic_field),
          twiss_(twiss),
          pipe_radius_(pipe_radius)
    {
    };
};

#endif // COOLER_H
