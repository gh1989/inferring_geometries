#include "pmcmc.h"
#include <getopt.h>
#include <stdio.h>

void output_time_series( Tensor<double, 2> chain, int trials, int parameters )
{
    FILE *fp;
    fp = fopen("output/pmcmc_timeseries.txt", "w+");    
    
    for( int i=0; i<trials; ++i )
    {
        for( int p=0; p<parameters; ++p )
            fprintf(fp, "%f\t", chain(i,p) );
        fprintf(fp, "\n");
    }

    fclose(fp);
}

int main(int argc, char *argv[])
{
    int trials = 1000;
    int rng_seed = 123456;
    int K = 25;  
    int N = 10;
    int T = 16;
    double dt = 0.001;
    double observation_noise_variance = 0.01;
    double trajectory_diffusion_sigma = 0.01;
    double proposal_c_variance = 0.05;
    const int parameters = 2;

    int opt;

    while( ( opt = getopt( argc, argv, ":T:R:K:N:L:d:o:t:c:" ) ) != EOF ) 
    {
    switch (opt)
        {
            case 'T':
            trials = atoi(optarg);
            break;
            case 'R':
            rng_seed = atoi(optarg);
            break;
            case 'K':
            K = atoi(optarg);
            break;
            case 'N':
            N = atoi(optarg);
            break;
            case 'L':
            T = atoi(optarg);
            break;
            case 'd':
            dt = atof(optarg);
            break;
            case 'o':
            observation_noise_variance = atof(optarg);
            break;  
            case 't':
            trajectory_diffusion_sigma = atof(optarg);
            break;  
            case 'c':
            proposal_c_variance = atof(optarg);
            break;      
        }
    }

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
    
    Tensor<double, 2> w(T, N);
    Tensor<double, 4> x(K, T, N, 2);
    Tensor<double, 3> y(K, T, 2);
    Tensor<double, 1> phat(trials);
    Tensor<double, 3> X(K, T, 2);
    Tensor<double, 3> X_star(K, T, 2);
    Tensor<double, 1> C(parameters);
    Tensor<double, 1> C_star(parameters);
    Tensor<double, 2> C_chain(trials, parameters);
    Tensor<double, 4> resampled( K, T, N, 2 );
    
    double acceptance=0;
    double log_acceptance_probability;
    double log_uniform_sample;
    double log_marginal_c;
    double log_marginal_c_star;
    double norm_total;
    double current_estimate_c;
    
    FourierSeries V(1);
    V.set_mode(1,1, 0.5 - 0.5*_Complex_I);
    
    generate_observations( r, y, 
                           K, T, V, dt, 
                           observation_noise_variance, 
                           trajectory_diffusion_sigma );

    for( int i=0; i<parameters; ++i )
    {
        C(i) = gsl_ran_gaussian( r, proposal_c_variance );
        C_chain(0,i) = C(i);
    }    

    V.set_mode(1, 1, C(0) + C(1)*_Complex_I);
    
    log_marginal_c = sequential_monte_carlo( r,
                                             x, w, y,
                                             phat,
                                             resampled,
                                             V,
                                             K, N, T, 
                                             dt,
                                             observation_noise_variance,
                                             trajectory_diffusion_sigma );

    double diff_total;
    for(int i=0; i<trials; ++i)
    {
        printf("%i...\n", i);
        
        for( int i=0; i<parameters; ++i )
            C_star(i) = gsl_ran_gaussian( r, proposal_c_variance ) + C(i);
        V.set_mode( 1,1, C_star(0) + C_star(1)*_Complex_I);

        log_marginal_c_star = sequential_monte_carlo( r,
                                                      x, w, y,
                                                      phat,
                                                      resampled,
                                                      V,
                                                      K, N, T, 
                                                      dt,
                                                      observation_noise_variance,
                                                      trajectory_diffusion_sigma );

        log_acceptance_probability = calculate_log_acceptance_probability( log_marginal_c,
                                                                           log_marginal_c_star, 
                                                                           C, 
                                                                           C_star, 
                                                                           proposal_c_variance,
                                                                           parameters );       
        
        log_uniform_sample = log( gsl_rng_uniform(r) );
        if ( log_uniform_sample < log_acceptance_probability )
        {         
            for( int ii=0; ii<parameters; ++ii ) C(ii) = C_star(ii);
            log_marginal_c = log_marginal_c_star;
            acceptance += 1;
            V.print_modes();
        }
        
        for( int ii=0; ii<parameters; ++ii ) C_chain(i, ii) = C(ii);
        
    }

    output_time_series( C_chain, trials, parameters );
    return 0;
}
