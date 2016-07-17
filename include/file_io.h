#pragma once
#include <complex.h>
#undef I
#include <math.h>
#include <stdio.h>
#include <Eigen/Dense>

using namespace Eigen;

void output_sigma_chain( int mcmc_trials, double *sigma_chain );
void output_average_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain );
void output_time_series_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain );
void output_posterior_density_complex( int M, int mcmc_trials, int param_size, double _Complex *chain, int bins);
void output_report( int M, int mcmc_trials, int param_size, double _Complex *chain, double _Complex *means, double acceptance_rate, 
                    int k, int P, double obs_sigma, double diff_sigma, double obs_delta ); // TODO: pass in settings/configuration object.