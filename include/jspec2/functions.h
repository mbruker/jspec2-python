#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <random>
#include <string>
#include <iostream>

bool iszero(double x);
bool iszero(double x, double err);
double rd(double x, double y, double z);
bool file_exists(std::string fileName);
std::string time_to_string();
std::string time_to_filename();

template <typename T>
void gaussian_random(int n, T& random_num, double sigma=1, double avg=0){
    std::default_random_engine generator;
    generator.seed(rand());
    std::normal_distribution<double> distribution(avg,sigma);
    for(int i=0; i<n; ++i) random_num[i] = distribution(generator);
}

template <typename T>
void gaussian_random_adjust(int n, T& random_num, double sigma, double avg=0) {
    double mean = 0;
    for(int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;

    double sigma_calc = 0;
    for(int i=0; i<n; ++i) {
        random_num[i] -= mean;
        sigma_calc += random_num[i]*random_num[i];
    }
    sigma_calc = sqrt(sigma_calc/n);

    double adjust_width = sigma/sigma_calc;
    for(int i=0; i<n; ++i) {
        random_num[i] *= adjust_width;
        random_num[i] += avg;
    }
}

template <typename T>
void uniform_random(int n, T& random_num, double r_min, double r_max){
    std::default_random_engine generator;
    generator.seed(rand());
    std::uniform_real_distribution<double> uniform_dis(r_min,r_max);
    for(int i=0; i<n; ++i) random_num[i] = uniform_dis(generator);
}

template <typename T>
void uniform_random_adjust(int n, T& random_num, double avg) {
    double mean = 0;
    for(unsigned int i=0; i<n; ++i) mean += random_num[i];
    mean /= n;
    double adjust = avg - mean;
    for(unsigned int i=0; i<n; ++i) random_num[i] += adjust;
}

template <typename T>
double rms(const T& v) {
    double sum = 0;
    for (auto _v : v) {
        sum += _v;
    }
    const double avg = sum/v.size();
    sum = 0;
    for (auto _v : v) {
        const double adj = _v - avg;
        sum += adj*adj;
    }
    return sqrt(sum/v.size());
}

#endif // FUNCTIONS_H
