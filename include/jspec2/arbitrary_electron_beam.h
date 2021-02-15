#ifndef ARBITRARY_ELECTRON_BEAM_H
#define ARBITRARY_ELECTRON_BEAM_H

#include <iostream>
#include <vector>

using std::vector;

typedef struct Box{
	double center[3] = {0,0,0};
	long int parent = 0;
	long int child[8] = {0};
	long int first_ptcl = 0;
	int n_ptcl = 0;
	int n_child = 0;
	double box_size = 0;
	long int first_ion = 0;
	int n_ion = 0;
} Box;

typedef struct Colleague{
	long int clg[28] = {0}; //At most 27 colleagues, the last 0 means the end.
} Colleague;

long int load_electrons(vector<double>& x, vector<double>& y, vector<double>& z, vector<double>& vx,
                                 vector<double>& vy, vector<double>& vz, std::string filename, long int n,
                                int skip = 0, bool binary = false, int n_buffer = 1000);

int create_e_tree(vector<double>& x, vector<double>& y, vector<double>& z, const unsigned long int n,
                  const unsigned int s, vector<Box> &tree, vector<unsigned long int>& list);

int create_ion_tree(double * x, double * y, double * z, const unsigned int n, vector<Box> &tree,
                vector<unsigned int>& list, unsigned int &idx_out);

std::ostream& operator<<(std::ostream& os, Box& box);
std::ostream& operator<<(std::ostream& os, Colleague& clg);

void create_e_tree(vector<double>& x, vector<double>& y, vector<double>& z, const long int n,
                  const int s, vector<Box> &tree, vector<long int>& list);

void create_ion_tree(vector<double>& x, vector<double>& y, vector<double>& z, const int n, vector<Box> &tree,
                vector<int>& list, int &idx_out);

void density(vector<Box> &tree, vector<long int>& list_e, vector<double>& vx, vector<double>& vy,
             vector<double>& vz, const long int ne,  vector<int>& list_i,
             int idx_out, const int ni, vector<double>& density_e,
             vector<double>& v_rms_t, vector<double>& v_rms_l) ;

void density(vector<Box> &tree, vector<long int>& list_e, vector<double>& vx, vector<double>& vy,
             vector<double>& vz, const long int ne,  vector<int>& list_i, int idx_out,
             const int ni, vector<double>& density_e, vector<double>& v_avg_z, vector<double>& v_rms_t,
             vector<double>& v_rms_l);


#endif // ARBITRARY_ELECTRON_BEAM_H
