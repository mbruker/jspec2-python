#ifndef TWISS_H
#define TWISS_H

#include <string>

class Twiss {
public:
    // Too lazy to write getters
    double bet_x = 0;
    double bet_y = 0;
    double alf_x = 0;
    double alf_y = 0;
    double disp_x = 0;
    double disp_y = 0;
    double disp_dx = 0;
    double disp_dy = 0;
    Twiss(double _bet_x,
          double _bet_y,
          double _alf_x,
          double _alf_y,
          double _disp_x,
          double _disp_y,
          double _disp_dx,
          double _disp_dy)
        : bet_x(_bet_x),
          bet_y(_bet_y),
          alf_x(_alf_x),
          alf_y(_alf_y),
          disp_x(_disp_x),
          disp_y(_disp_y),
          disp_dx(_disp_dx),
          disp_dy(_disp_dy)
    {
        // TODO error handling
    }
    std::string toString() const
    {
        return "beta=(" + std::to_string(bet_x) + " ; " + std::to_string(bet_y) + ")\n"
             + "alpha=(" + std::to_string(alf_x) + " ; " + std::to_string(alf_y) + ")\n"
             + "disp=(" + std::to_string(disp_x) + " ; " + std::to_string(disp_y) + ")\n"
             + "disp_d=(" + std::to_string(disp_dx) + " ; " + std::to_string(disp_dy) + ")";
    }
};

#endif
