#pragma once
#include <complex.h>
#undef I
#include <math.h>
#include <stdio.h>
#include <Eigen/Dense>

using namespace Eigen;

void output_average_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain );
void output_time_series_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain );
void output_posterior_density_complex( int M, int mcmc_trials, int param_size, double _Complex *chain, int bins);