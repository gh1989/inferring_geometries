#include "pmcmc.h"

int main()
{
    // MCMC trials
    const int trials = 10000;
    // Random number generator seed
    const int rng_seed = 123;
    // Number of particles
    const int N = 1;
    // Time horizon
    const int T = 16;
    // Diffusion time-step
    const double dt = 0.001;
    // Observation noise variance, y_t = x_t + e_t, with e_t ~ N(0, var)
    const double observation_noise_variance = 0.1;
    // The Langevin diffusion matrix coefficient
    const double trajectory_diffusion_sigma = 0.001;
    // The proposal for the parameter variance for Gaussian
    const double proposal_c_variance = 0.05;
    // The number of parameters
    const int parameters = 1;
  
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
    
    Tensor<double, 2> w(T, N);
    Tensor<double, 3> x(T, N, 2);
    Tensor<double, 2> y(T, 2);
    Tensor<double, 1> phat(trials);
    Tensor<double, 2> X(T, 2);
    Tensor<double, 2> X_star(T, 2);
    Tensor<double, 1> C(parameters);
    Tensor<double, 1> C_star(parameters);
    Tensor<double, 2> C_chain(trials, parameters);
    
    double log_acceptance_probability;
    double log_uniform_sample;
    double marginal_likelihood_c;
    double marginal_likelihood_c_star;
    double real_c = 0.5;
    
    FourierSeries V(1);
    V.set_mode(1,1,real_c);
    
    generate_observations( r, y, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    C(0) = gsl_ran_gaussian( r, proposal_c_variance );
    C_chain(0,0) = C(0);
    V.set_mode(1,1,C(0));
    marginal_likelihood_c = sequential_monte_carlo( r, x, w, y, phat, V, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
    sample_smc_approximation( r, N, T, x, w, X );
    
    for(int i=0; i<trials; ++i)
    {
        C_star(0) = gsl_ran_gaussian( r, proposal_c_variance ) + C(0);
        V.set_mode(1, 1, C_star(0));
        marginal_likelihood_c_star = sequential_monte_carlo( r, x, w, y, phat, V, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
        sample_smc_approximation( r, N, T, x, w, X_star);

        log_acceptance_probability = calculate_log_acceptance_probability( marginal_likelihood_c, marginal_likelihood_c_star, C, C_star, proposal_c_variance);
        log_uniform_sample = log( gsl_rng_uniform(r) );
        if ( log_uniform_sample < log_acceptance_probability )
        {
            C(0) = C_star(0);
        }
        C_chain(i) = C(0);
        printf("C_chain(%i) = %f. \n", i, C_chain(i) );
    }
    
    return 0;
}