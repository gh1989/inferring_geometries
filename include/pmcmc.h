#pragma once

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

#include "smc.h"
#include "fourier_series.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void sample_smc_approximation( gsl_rng *r,
                               int K, 
                               int N,
                               int T,
                               Tensor<double, 4> &x, 
                               Tensor<double, 2> &w,
                               Tensor<double, 3> &X );
                               
double calculate_log_acceptance_probability( double log_marginal_likelihood_c,
                                             double log_marginal_likelihood_c_star,
                                             Tensor<double, 1> &C, 
                                             Tensor<double, 1> &C_star,
                                             double proposal_c_variance );
