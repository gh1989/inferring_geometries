#include "file_io.h"
#include "mcmc.h"

int main( int argc, char *argv[] )
{
	// s: The random number generator seed for GSL library.
	int rng_seed = 123456789;
	// P: Number of timesteps per path generated for the particle.
	int path_length = 64;
	// N: Iterations in the MCMC algo.
	int mcmc_trials = 1000;
	// c: Parameter perturbation amount.
	double param_delta = 0.05;
	// d: Observation Y perturbation amount.
	double obs_sigma = 0.01;
	// j: observational timestep
	double obs_delta = 0.008;
	// p: The timestep for the particle's path in the discretisation.
	double path_delta = 0.001;
	// b: Burn in time.
	int burn = 0;
	// o: sigma of the diffusion
	double real_sigma = 0.01;
	// k: number of parallel particles
	int parallel_paths = 20;
	// h: histogram bins
	int bins = 100;
	// x: sigma being a parameter to be inferred the lognormal proposal needs a variance
	double param_delta_sigma = 0.05;
	// z: zeta for the lognormal proposal
	double zeta = log( 0.01 );
	// X: infer diffusion coefficient?
	bool infer_sigma = false;	
	// The option.
	int opt;
	
	// Get options from command line.
	while( ( opt = getopt( argc, argv, ":s:P:N:c:d:p:b:q:o:k:h:j:x:z:X:" ) ) != EOF ) 
	{
	switch (opt)
		{
			case 's':
			rng_seed = atoi(optarg);
			break;
			
			case 'P':
			path_length = atoi(optarg);
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
			real_sigma = atof(optarg);
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

			case 'x':
			param_delta_sigma = atof(optarg);
			break;
			
			case 'z':
			zeta = atof(optarg);
			break;
			
			case 'X':
			infer_sigma = atoi(optarg);
			break;
		}
	}
	
	// Expansion cutoff and implied param_size ( including (0,0) -> 0 ).
	int M = 1;
	FourierSeries v(M);
	int acceptance_frequency = 0;
	double alpha;
	double acceptance_rate;
	double _Complex average_constant = 0;
	double p_param_prob,g_param_prob,p_paths_prob;
	double log_u, log_alpha;
	
	printf("Path delta set to: %f. \n", path_delta );
	printf("Observation delta set to: %f. \n", obs_sigma );
	
	FourierSeries real_v(M);
	real_v.set_mode( 1, 1, 0.5 - 0.5*_Complex_I);
		
	// Identity matrix for the diffusion matrix.
	Matrix2d C;
	C << real_sigma, 0, 0, real_sigma;
	
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
	Vector2d *real_path =  (Vector2d*) malloc( parallel_paths*path_length*sizeof( Vector2d ) ); 
	Vector2d *x_star = (Vector2d*) malloc( path_length*parallel_paths*sizeof( Vector2d ) );
	
	int observed_length;
	int actual_per_observed;
	
	observed_length = path_length * path_delta / obs_delta;
	actual_per_observed = obs_delta / path_delta;
	
	printf("[Debug] observed_length: %i, actual_per_observed: %i. \n", observed_length, actual_per_observed );
	
	Vector2d *x = (Vector2d*) malloc( path_length*parallel_paths*sizeof( Vector2d ) );
	Vector2d *y = (Vector2d*) malloc( observed_length*parallel_paths*sizeof( Vector2d ) );
	
	double _Complex c;
	double _Complex c_star; 
	double _Complex *chain = (double _Complex*) malloc(mcmc_trials*sizeof(double _Complex));
	
	
	for( int i = 0; i<parallel_paths; ++i)
		real_starting_points[i] << gsl_rng_uniform(r)*2-1, gsl_rng_uniform(r)*2-1;
		
	for( int i=0; i<parallel_paths; ++i )
		em_overdamped_langevin( path_length, path_delta, C, &real_v, real_path+i*path_length, r, real_starting_points[i]);

	for( int i=0; i<parallel_paths; ++i)
	{
		for( int j=0; j<observed_length; ++j )
		{
			noise << gsl_ran_gaussian( r, obs_sigma ), gsl_ran_gaussian( r, obs_sigma );
			y[i*observed_length+j] = real_path[i*path_length+j*actual_per_observed] + noise;
		}
		observed_starting_points[i] = y[i*observed_length];
	}
	
	double *sigma_chain; 
	double sigma;
	double sigma_star;
	if( infer_sigma )
	{
		sigma_chain = (double*)malloc( sizeof(double)*mcmc_trials );
		sigma = gsl_ran_lognormal(r, zeta, param_delta_sigma );
		C << sigma, 0, 0, sigma;
		sigma_chain[0] = sigma;
	}
	c = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I;
	chain[0] = c;
	v.set_mode( 1, 1, c );

	for( int i=0; i<parallel_paths; ++i )
		em_overdamped_langevin( path_length, path_delta, C, &v, x+i*path_length, r, observed_starting_points[i]);
	
	printf("\n ********** \n Starting MCMC. \n ********** \n");
	for( int i=0; i<mcmc_trials+burn; ++i )
	{
		c_star = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I + c; // c* ~ N(c, param_delta).
		
		if( infer_sigma )
		{
			sigma_star = gsl_ran_lognormal(r, log(sigma), param_delta_sigma );
			C << sigma_star, 0, 0, sigma_star;
		}
		v.set_mode( 1, 1, c_star );
		for( int ii=0; ii<parallel_paths; ++ii )
			em_overdamped_langevin( path_length, path_delta, C, &v, x_star+ii*path_length, r, observed_starting_points[ii]); // Generate x*
		
		log_u = log( gsl_rng_uniform(r) );
		p_param_prob = log_p( c_star, param_delta ) - log_p( c, param_delta );		
		
		if( infer_sigma )
		{
			p_param_prob +=  log_p_sigma( sigma_star, param_delta_sigma ) -  log_p_sigma( sigma, param_delta_sigma );
		}
		
		//g_param_prob = log_g( c, c_star, param_delta, param_size) - log_g( c_star, c, param_delta, param_size);
		p_paths_prob = log_p( y, x_star, observed_length*parallel_paths, obs_sigma, actual_per_observed) - log_p(y, x, observed_length*parallel_paths, obs_sigma, actual_per_observed);	
		log_alpha = p_param_prob + g_param_prob + p_paths_prob;
		
		if (log_u<log_alpha)
		{
			//printf("[MCMC] n=%i.\n -> Accepted: c_star! \n", i);
			if (i>burn) acceptance_frequency += 1;
			c = c_star;
			sigma = sigma_star;
			for(int j=0; j<path_length*parallel_paths; ++j) x[j] = x_star[j];
			//v.print_modes();
		}
		 if (i>burn) {
			chain[ (i-burn) ] = c;
			if( infer_sigma ) sigma_chain[ (i-burn) ] = sigma;
		}
	}
	acceptance_rate = double( acceptance_frequency ) / mcmc_trials;
	double _Complex averages[1];
	
	printf("\n ********** \n Chain completed. \n ********** \n");
	printf("Acceptance rate: %f \n", acceptance_rate );
	
	int ii, jj, idx;
	printf("[MCMC: Terminal constants] ");
	v.set_mode( 1, 1, c );
	v.print_modes();
	
	for( int i=0; i<mcmc_trials; ++i )
		averages[0] += chain[i];
	averages[0] /= mcmc_trials;
		
	printf("[MCMC: Averages]  " );		
	v.set_mode( 1,1, averages[0] );
	v.print_modes();
	
	if( infer_sigma ) output_sigma_chain( mcmc_trials, sigma_chain ); 
	output_average_file_complex(M, mcmc_trials, 1, chain);
	output_time_series_file_complex(M, mcmc_trials, 1, chain);
	output_posterior_density_complex(M, mcmc_trials, 1, chain, bins);
	output_report(  M, mcmc_trials, 1, chain, averages, acceptance_rate, parallel_paths, path_length, obs_sigma, sigma, obs_delta );
	//free(real_path);
	//free(x);
	//free(x_star);	
	//free(chain);
}