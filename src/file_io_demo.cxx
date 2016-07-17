#include "file_io.h"
#include <complex.h>
#undef I
#include <math.h>
#include <stdio.h>

int main()
{
    const int M = 1;
    const int mcmc_trials = 16;
    const int param_size = 1;
    const int bins = 32;
    
    double acceptance_rate = 0.52;
    double _Complex chain[mcmc_trials*param_size];
    double _Complex means[param_size];
    
    for( int j = 0; j < param_size; ++j )
    {
        means[j] = 0.25;
        for( int i = 0; i < mcmc_trials; ++i )
        {
            chain[ i*param_size + j ] = i;
        }
    }
    int k = 20;
    int P = 16;
    double obs_sigma = 0.15; 
    double diff_sigma = 0.1;
    double obs_delta = 0.1;
    
    output_average_file_complex( M, mcmc_trials, param_size, chain );
    output_time_series_file_complex( M, mcmc_trials, param_size, chain );
    output_posterior_density_complex( M, mcmc_trials, param_size, chain, bins);
    output_report( M, mcmc_trials, param_size, chain, means, acceptance_rate, k, P, obs_sigma, diff_sigma, obs_delta );

    return 0;
}