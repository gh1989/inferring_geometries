#pragma once
#include <Eigen/Dense>
#include "euler_maruyama.h"
#include "fourier_series.h"
#include <gsl/gsl_cdf.h>
#include <math.h>
#include <stdio.h>

double log_p_sigma( double sigma, double param_delta );
double log_p( double c, double param_delta);
double log_p( double _Complex c, double param_delta );
double log_p( double _Complex *c, double param_delta, int length );
double log_g( double c_star, double c, double param_delta );
double log_g_sigma( double sigma_star, double sigma, double param_delta );
double log_g( double _Complex *c_star, double _Complex *c, double param_delta, int length );
double log_p( Vector2d *y, Vector2d *x, int path_steps, double obs_delta, int actual_per_observed );