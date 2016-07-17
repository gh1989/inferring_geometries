#pragma once
#include <assert.h>
#include <fftw3.h>
#include "fourier_series.h"

struct Grid
{
    Grid( double, double, double, double, int, int );
    ~Grid();

    double min_x;
    double min_y;
    double max_x;
    double max_y;
    double x_delta;
    double y_delta;
    
    double Lx;
    double Ly;
    
    int size_x;
    int size_y;
    
    int Nx;
    int Ny;
    
    Array<double, Dynamic, 1> discrete_domain_x;
    Array<double, Dynamic, 1> discrete_domain_y;
    Array<double, Dynamic, Dynamic> discrete_range;
    Array<double, Dynamic, Dynamic> grad_discrete_range;
    
    Vector2d interpolate_gradient( double x, double y );
    Vector2d interpolate_gradient( Vector2d v );
    double evaluate( double x, double y );
    double evaluate( Vector2d v );
    void discretise( FourierSeries *pf_series );

    fftw_complex *y;
    fftw_complex *Y;
    fftw_complex *vx;
    fftw_complex *vy;
    fftw_complex *Vx;
    fftw_complex *Vy;
    
};