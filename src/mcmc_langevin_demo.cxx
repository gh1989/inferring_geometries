#include "file_io.h"
#include "mcmc.h"
#include <getopt.h>

int main( int argc, char *argv[] )
{
    // s: The random number generator seed for GSL library.
    int rng_seed = 123456789;
    // P: Number of timesteps per path generated for the particle.
    int path_steps = 16;
    // N: Iterations in the MCMC algo.
    int mcmc_trials = 1000;
    // c: Parameter perturbation amount.
    double param_delta = 0.05;
    // d: Observation Y perturbation amount.
    double obs_sigma = 0.01;
    // j: observational timestep
    double obs_delta = 0.001;
    // p: The timestep for the particle's path in the discretisation.
    double path_delta = 0.001;
    // b: Burn in time.
    int burn = 0;
    // o: sigma of the diffusion
    double sigma = 0.01;
    // k: number of parallel particles
    int parallel_paths = 20;
    // h: histogram bins
    int bins = 100;
    // The option.
    int opt;
    
    // Get options from command line.
    while( ( opt = getopt( argc, argv, ":s:P:N:c:d:p:b:q:o:k:h:j:" ) ) != EOF ) 
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
            obs_sigma = atof(optarg);
            break;
            
            case 'p':
            path_delta = atof(optarg);
            break;
            
            case 'b':
            burn = atoi(optarg);
            break;
            
            case 'o':
            sigma = atof(optarg);
            break;
            
            case 'k':
            parallel_paths = atoi(optarg);
            break;
            
            case 'h':
            bins = atoi(optarg);
            break;
            
            case 'j':
            obs_delta = atof(optarg);
            break;

        }
    }
    
    // Expansion cutoff and implied param_size ( including (0,0) -> 0 ).
    int M = 1;
    int param_size = 2*M*(M+1);
    FourierSeries v(1);
    int acceptance_frequency = 0;
    double alpha;
    double acceptance_rate;
    double _Complex average_constant = 0;
    double p_param_prob,g_param_prob,p_paths_prob;
    double log_u, log_alpha;
    
    printf("Path delta set to: %f. \n", path_delta );
    printf("Observation delta set to: %f. \n", obs_sigma );
    
    FourierSeries real_v(M);
    real_v.set_mode( 0, 1, 0.5 - 0.5*_Complex_I);
    real_v.set_mode( 1, 0, 0.5 - 0.5*_Complex_I);
        
    // Identity matrix for the diffusion matrix.
    Matrix2d C;
    C << sigma, 0, 
         0, sigma;
    
    // Setup the random number generator
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
    
    // Some useful temporary vectors for random starting point, observational noise...
    Vector2d noise;
    Vector2d random_starting_point;
    
    // The real starting points of the parallel observations, we do not want the particles to be stuck in "modes" of the
    // potential V so this helps.
    Vector2d *real_starting_points = (Vector2d*) malloc( parallel_paths*sizeof( Vector2d ) );    
    Vector2d *observed_starting_points = (Vector2d*) malloc( parallel_paths*sizeof( Vector2d ) );
    Vector2d *real_path =  ( Vector2d* ) malloc( parallel_paths*path_steps*sizeof( Vector2d ) ); 
    Vector2d *x_star = (Vector2d*) malloc( path_steps*parallel_paths*sizeof( Vector2d ) );
    Vector2d *x = (Vector2d*) malloc( path_steps*parallel_paths*sizeof( Vector2d ) );
    Vector2d *y = ( Vector2d* ) malloc( parallel_paths*path_steps*sizeof( Vector2d ) );
    
    double _Complex *c = (double _Complex*) malloc(param_size*sizeof(double _Complex));
    double _Complex *c_star = (double _Complex*) malloc(param_size*sizeof(double _Complex));
    double _Complex *chain = (double _Complex*) malloc(mcmc_trials*param_size*sizeof(double _Complex));
    
    for( int i = 0; i<parallel_paths; ++i)
        real_starting_points[i] << gsl_rng_uniform(r)*2-1, gsl_rng_uniform(r)*2-1;
        
    for( int i=0; i<parallel_paths; ++i )
        em_overdamped_langevin( path_steps, path_delta, C, &real_v, real_path+i*path_steps, r, real_starting_points[i]);

    for( int i=0; i<parallel_paths; ++i)
    {
        for( int j=0; j<path_steps; ++j )
        {
            noise << gsl_ran_gaussian( r, obs_sigma ), gsl_ran_gaussian( r, obs_sigma );
            y[i*path_steps+j] = real_path[i*path_steps+j] + noise;
        }
        observed_starting_points[ i ] = y[i*path_steps];
    }
    
    for( int j=0; j<param_size; j++ )
    {
        c[j] = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I;
        chain[j] = c[j];
    }
    
    v.set_modes( c );
    
    for( int i=0; i<parallel_paths; ++i )
        em_overdamped_langevin( path_steps, path_delta, C, &v, x+i*path_steps, r, observed_starting_points[i]);
    
    printf("\n ********** \n Starting MCMC. \n ********** \n");
    for( int i=0; i<mcmc_trials+burn; ++i )
    {
        // Propose constants
        for( int j=0; j<param_size; ++j)
            c_star[j] = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I + c[j]; // c* ~ N(c, param_delta).
        
        v.set_modes( c_star );
        for( int ii=0; ii<parallel_paths; ++ii )
            em_overdamped_langevin( path_steps, path_delta, C, &v, x_star+ii*path_steps, r, observed_starting_points[ii]); // Generate x*
        
        log_u = log( gsl_rng_uniform(r) );
        p_param_prob = log_p( c_star, param_delta, param_size) - log_p( c, param_delta, param_size);        
        //g_param_prob = log_g( c, c_star, param_delta, param_size) - log_g( c_star, c, param_delta, param_size);
        p_paths_prob = log_p( y, x_star, path_steps*parallel_paths, obs_sigma, 1) - log_p(y, x, path_steps*parallel_paths, obs_sigma, 1);    
        log_alpha = p_param_prob + g_param_prob + p_paths_prob;
        
        if (log_u<log_alpha)
        {
            //printf("[MCMC] n=%i.\n -> Accepted: c_star! \n", i);
            if (i>burn) acceptance_frequency += 1;
            for( int j=0; j<param_size; ++j) c[j] = c_star[j];
            for(int j=0; j<path_steps*parallel_paths; ++j) x[j] = x_star[j];
            //v.print_modes();
        }
         if (i>burn) {
            for( int j=0; j<param_size; ++j)
                chain[ (i-burn)*param_size + j ] = c[j];
        }
    }
    acceptance_rate = double( acceptance_frequency ) / mcmc_trials;
    double _Complex averages[param_size];
    
    printf("\n ********** \n Chain completed. \n ********** \n");
    printf("Acceptance rate: %f \n", acceptance_rate );
    
    int ii, jj, idx;
    printf("[MCMC: Terminal constants] ");
    v.set_modes( c );
    v.print_modes();
    
    for( int i=0; i<mcmc_trials; ++i )
        for( int j=0; j<param_size; ++j )
            averages[j] += chain[i*param_size + j];
    for( int k=0; k<param_size; ++k )
        averages[k] /= mcmc_trials;
        
    printf("[MCMC: Averages]  " );        
    v.set_modes( averages );
    v.print_modes();
    
    output_average_file_complex(M, mcmc_trials, param_size, chain);
    output_time_series_file_complex(M, mcmc_trials, param_size, chain);
    output_posterior_density_complex(M, mcmc_trials, param_size, chain, bins);
    //free(real_path);
    //free(x);
    //free(x_star);    
    //free(chain);
    
}
