#ifndef RING_H
#define RING_H

#include <cassert>
#include <string>
#include <vector>

using std::vector;

class IonBeam;

class Lattice{
    vector<double> s_;
    vector<double> betx_;
    vector<double> alfx_;
    vector<double> mux_;
    vector<double> dx_;
    vector<double> dpx_;
    vector<double> bety_;
    vector<double> alfy_;
    vector<double> muy_;
    vector<double> dy_;
    vector<double> dpy_;
    vector<double> l_element_;
    int n_element_;
    double circ_;
 public:
    double s(int i) const {return s_.at(i);}
    double betx(int i) const {return betx_.at(i);}
    double alfx(int i) const {return alfx_.at(i);}
    double mux(int i) const {return mux_.at(i);}
    double dx(int i) const {return dx_.at(i);}
    double dpx(int i) const {return dpx_.at(i);}
    double bety(int i) const {return bety_.at(i);}
    double alfy(int i) const {return alfy_.at(i);}
    double muy(int i) const {return muy_.at(i);}
    double dy(int i) const {return dy_.at(i);}
    double dpy(int i) const {return dpy_.at(i);}
    int n_element() const {return n_element_;}
    double l_element(int i) const {return l_element_.at(i);}
    double circ() const {return circ_;}
    Lattice(std::string filename);
    Lattice(double circ);
};

class Ring {
protected:
    double beta_s_ = 0;         //Synchrotron function, use to calculate rms bunch length from momentum spread
    double f0_ = 0;               // revolution frequency.
    double w0_ = 0;       // angular revolution frequency.
    const IonBeam *ionBeam_;
    Lattice lattice_;
    // Tunes
    double qx_ = 0;
    double qy_ = 0;
    double qs_ = 0;
    // RF
    double rf_voltage_ = 0;     // RF voltage in V
    int rf_h_ = 1;               // RF harmonic number
    double rf_phi_ = 0;          // RF phase in 2*pi
    double gamma_tr_ = 0;        // Transition gamma
 public:
    double beta_s() const { return beta_s_; }
    double circ() const { return lattice_.circ(); }
    double f0() const {return f0_;}
    double w0() const {return w0_;}
    double rf_voltage() const { return rf_voltage_; }
    double rf_h() const { return rf_h_; }
    double rf_phi() const { return rf_phi_; }
    double gamma_tr() const { return gamma_tr_; }
    double slip_factor() const;
//    double calc_sync_tune_by_rf() const;
    const IonBeam& ionBeam() const {return *ionBeam_;}
    const Lattice& lattice() const {return lattice_;}
    double qx() const { return qx_; }
    double qy() const { return qy_; }
    double qs() const { return qs_; }
    Ring(const Lattice &lattice,
         const IonBeam *beam,
         double qx = 0,
         double qy = 0,
         double qs = 0,
         double rf_voltage = 0,
         int rf_h = 1,
         double rf_phi = 0,
         double gamma_tr = 0
        );
    void update_bet_s();
    void update_rf_voltage();
};

#endif // RING_H
