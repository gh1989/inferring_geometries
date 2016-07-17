#include "pmcmc.h"

int main()
{
    // MCMC trials
    const int trials = 10000;
    // Random number generator seed
    const int rng_seed = 123;
    // Number of particles
    const int N = 8;
    // Time horizon
    const int T = 100;
    // Diffusion time-step
    const double dt = 0.1;
    // Observation noise variance, y_t = x_t + e_t, with e_t ~ N(0, var)
    const double observation_noise_variance = 0.1;
    // The Langevin diffusion matrix coefficient
    const double trajectory_diffusion_sigma = 0.1;
  
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
    
    Tensor<double, 2> w(T, N);
    Tensor<double, 3> x(T, N, 2);
	Tensor<double, 2> y(T, 2);
	    
    FourierSeries V(1);
    V.set_mode(1,1,0.5);
        
    sequential_monte_carlo( r, x, w, y, V, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
    return 0;
}