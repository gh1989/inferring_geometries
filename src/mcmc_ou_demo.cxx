#include "mcmc.h"

int main( int argc, char *argv[] )
{
    // s: The random number generator seed for GSL library.
    int rng_seed = 14081989;
    // P: Number of timesteps per path generated for the particle.
    int path_steps = 10;
    // N: Iterations in the MCMC algo.
    int mcmc_trials = 10;
    // c: Parameter perturbation amount.
    double param_delta = 1.0;
    // d: Observation Y perturbation amount.
    double obs_delta = 0.1;
    // p: The timestep for the particle's path in the discretisation.
    double path_delta = 0.001;
    // b: Burn in time.
    int burn = 0;
    // q: The real parameter to be inferred.
    double real_parameter_c = 0.5;

    // The option.
    int opt;
    
    // Get options from command line.
    while( ( opt = getopt( argc, argv, ":s:P:N:c:d:p:b:q:" ) ) != EOF ) 
    {
    switch (opt)
        {
            case 's':
            rng_seed = atoi(optarg);
            break;
            
            case 'P':
            path_steps = atoi(optarg);
            break;
            
            case 'N':
            mcmc_trials = atoi(optarg);
            break;
            
            case 'c':
            param_delta = atof(optarg);
            break;
            
            case 'd':
            obs_delta = atof(optarg);
            break;
            
            case 'p':
            path_delta = atof(optarg);
            break;
            
            case 'b':
            burn = atoi(optarg);
            break;
            
            case 'q':
            real_parameter_c = atof(optarg);
            break;
        }
    }
    
    printf("Constant set to: %f. \n", real_parameter_c );
    printf("Path delta set to: %f. \n", path_delta );
    printf("Observation delta set to: %f. \n", obs_delta );
    
    // Setup the random number generator
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
    
    // Generate a path and the observation.
    Vector2d *real_path, *y;
    real_path = ( Vector2d* ) malloc( path_steps*sizeof( Vector2d ) );
    y = ( Vector2d* ) malloc( path_steps*sizeof( Vector2d ) );
    em_ornstein_uhlenbeck( path_steps, path_delta, real_parameter_c, real_path, r );
    
    Vector2d noise;
    for( int i=0; i<path_steps; ++i )
    {
        noise << gsl_ran_gaussian( r, obs_delta ), 
                 gsl_ran_gaussian( r, obs_delta );
        
        y[i] = real_path[i] + noise;
    }

    // Attempt to obtain the value of real_parameter_c with mcmc algo.
    int acceptance_frequency = 0;
    double alpha;
    double acceptance_rate;
    double average_constant = 0;
    double c, c_star;
    double *chain = (double*) malloc(mcmc_trials*sizeof(double));
    double u, log_alpha;
    
    c = gsl_ran_gaussian( r, param_delta );
    chain[0] = c;
    Vector2d *x_star = (Vector2d*) malloc( path_steps*sizeof( Vector2d ) );
    Vector2d *x = (Vector2d*) malloc( path_steps*sizeof( Vector2d ) );
    
    em_ornstein_uhlenbeck( path_steps, path_delta, c, x, r );
    
    for( int i=0; i<mcmc_trials+burn; ++i )
    {
        c_star = gsl_ran_gaussian( r, param_delta ) + c; // c* ~ N(c, param_delta).
        em_ornstein_uhlenbeck( path_steps, path_delta, c_star, x_star, r); // Generate x*
        
        u = gsl_rng_uniform(r);
        log_alpha = log_p( c_star, param_delta) - log_p( c, param_delta);
        //log_alpha += log_g( c, c_star, param_delta) - log_g( c_star, c, param_delta);
        log_alpha += log_p( y, x_star, path_steps, obs_delta, 1) - log_p(y, x, path_steps, obs_delta, 1);
        alpha = exp( log_alpha );
        if (u<alpha)
        {
            // Accept (c*, x*)
            if (i>burn)    {
            acceptance_frequency += 1;
            chain[ i+1-burn ] = c_star;
            }
            c = c_star;
            
            for(int j=0; j<path_steps; ++j)
                x[j] = x_star[j];
        }
        else
        {
            // Reject (c*, x*)
            if (i>burn) chain[ i+1-burn ] = c;
        }
    }
    acceptance_rate = double( acceptance_frequency ) / mcmc_trials;
    
    printf("Chain completed...");
    printf("Acceptance rate: %f \n", acceptance_rate );
    printf("Parameter ended with: %f \n", c_star );
    
    for( int i=0; i<mcmc_trials; ++i )
        average_constant += chain[i];
    average_constant /= mcmc_trials;
    
    printf("Average value for constant from chain: %f", average_constant);
    
    // File IO
    int data_points = 1000;
    int block_size = mcmc_trials / data_points;
    int left_over_size = mcmc_trials % data_points;
    
    double cumulative_sum = 0.0;
    
    // Moving average data.
    FILE *fp;
    fp = fopen("Testing/MCMC/data.txt", "w+");
    fprintf(fp, "#n\tmean\n");
    
    for(int i=0; i<data_points; ++i)
    {
        for(int j=i*block_size; j<(i+1)*block_size; ++j)
            cumulative_sum += chain[j];
        fprintf(fp, "%i\t%f\n", (i+1)*block_size, cumulative_sum/(i+1)/block_size );
    }

    for(int i=mcmc_trials-left_over_size; i<mcmc_trials; i++)
        cumulative_sum += chain[i];    
    fprintf(fp, "%i\t%f\n", mcmc_trials, cumulative_sum/mcmc_trials );
    fclose(fp);
    
    // Time series data.
    fp = fopen("Testing/MCMC/data_ts.txt", "w+");
    fprintf(fp, "#n\tvalue\n");
    
    for(int i=0; i<mcmc_trials; ++i)
        fprintf(fp, "%i\t%f\n", i, chain[i] );
    
    fclose(fp);
        
    // free allocated memory
    free(real_path);
    free(x);
    free(x_star);    
    free(chain);
}