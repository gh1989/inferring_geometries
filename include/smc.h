#pragma once

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

#include "fourier_series.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace Eigen;

void generate_observations( gsl_rng *r, 
                            Tensor<double, 3> &y, 
                            int K, int T, FourierSeries &V, double dt,
                            double observation_noise_variance, double trajectory_diffusion_sigma  );

void generate_particle_samples( gsl_rng *r, 
                                Tensor<double, 4> &x, 
                                Tensor<double, 3> &y,
                                int K, int N, int t, FourierSeries &V, double dt, 
                                double observation_noise_variance, double trajectory_diffusion_sigma );

void generate_particle_samples( gsl_rng *r, 
                                Tensor<double, 4> &x, 
                                Tensor<double, 3> &y,
                                int K, int N, int t, int M,
                                FourierSeries &V, 
                                double dt, double ds,
                                double observation_noise_variance, 
                                double trajectory_diffusion_sigma );
                                
void simulate_forwards( gsl_rng*, Tensor<double, 3>&, FourierSeries& );

double assign_weights( Tensor<double, 4> &x, 
                     Tensor<double, 2> &w,
                     Tensor<double, 3> &y,
                     int K, int N, int t,
                     double observation_noise_variance );

double sequential_monte_carlo( gsl_rng *r,
                               Tensor<double, 4> &x,
                               Tensor<double, 2> &w,
                               Tensor<double, 3> &y,
                               Tensor<double, 1> &phat,
                               Tensor<double, 4> &resampled,
                               FourierSeries &V,
                               int K, int N, int T,
                               double dt, 
                               double observation_noise_variance,
                               double trajectory_diffusion_sigma );
    
// Diffusion Bridge SMC    
double sequential_monte_carlo( gsl_rng *r,
                               Tensor<double, 4> &x,
                               Tensor<double, 2> &w,
                               Tensor<double, 3> &y,
                               Tensor<double, 1> &phat,
                               Tensor<double, 4> &resampled,
                               FourierSeries &V,
                               int K, int N, int T, int M,
                               double dt, double ds,
                               double observation_noise_variance,
                               double trajectory_diffusion_sigma );
                      
double estimate_marginal_likelihood( int t, int N, 
                                     Tensor<double, 2> &w, 
                                     Tensor<double, 1> &phat );
                                     
void resample( gsl_rng *r,
               gsl_ran_discrete_t *g,
               Tensor<double,4> &x,
               Tensor<double,2> &w,
               Tensor<double, 4> &resampled,
               int K, int N, int T );

Vector2d bridge_drift(  int j, 
                        double ds,
                        Vector2d xj, 
                        Vector2d y_end,
                        FourierSeries &V,
                        double observation_noise_variance, 
                        double trajectory_diffusion_sigma);
                        
double bridge_variance( int j,
                        int ds,
                        double observation_noise_variance, 
                        double trajectory_diffusion_sigma );