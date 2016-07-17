#pragma once

#include <Eigen/Dense>
#include <complex.h>
#undef I

using namespace Eigen;

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

// M is the cutoff for |i| <= M.
// A full set -M <= i = M is 2*M+1 values
// Due to the reality constraint v*_{-k} = v_{k}.
const double _Complex two_pi_i = 2*M_PI*_Complex_I;

struct FourierSeries
{
    Vector2d grad( Vector2d );
    Vector2d grad( double, double );
    double evaluate( Vector2d );
    double evaluate( double, double );
    void set_mode( int,  int, double _Complex );
    void set_modes( double _Complex* );
    double _Complex get_mode( int, int );
    void get_modes( double _Complex* );
    FourierSeries(int M_);
    ~FourierSeries();
    double _Complex *modes;
    int total_modes;
    void print_modes();
    int M;
};