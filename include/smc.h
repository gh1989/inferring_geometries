#pragma once
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>
#include "fourier_series.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Eigen;

void generate_observations( gsl_rng *r, 
							Tensor<double, 2> &y, 
							int T, FourierSeries &V, double dt,
							double observation_noise_variance, double trajectory_diffusion_sigma  );

void generate_particle_samples( gsl_rng*, Tensor<double, 3>&, int, int );


void assign_weights( Tensor<double, 3> &x, 
					 Tensor<double, 2> &w,
					 Tensor<double, 2> &y,
					 int N, int t,
	                 double observation_noise_variance );

double calculate_weight( Tensor<double,2> &y, Tensor<double,3> &x, 
						 int N, int i, int t,
						 double observation_noise_variance );

