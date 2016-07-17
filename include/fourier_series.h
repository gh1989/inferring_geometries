#pragma once

#include <complex.h>
#undef I

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>

using namespace Eigen;

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