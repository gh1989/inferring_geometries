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
	double obs_delta = 0.01;
	// p: The timestep for the particle's path in the discretisation.
	double path_delta = 0.001;
	// b: Burn in time.
	int burn = 0;
	// o: sigma of the diffusion
	double sigma = 0.1;
	// k: number of parallel particles
	int parallel_paths = 1;
	
	// The option.
	int opt;
	
	// Get options from command line.
	while( ( opt = getopt( argc, argv, ":s:P:N:c:d:p:b:q:o:k:" ) ) != EOF ) 
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
			
			case 'o':
			sigma = atof(optarg);
			break;
			
			case 'k':
			parallel_paths = atoi(optarg);
			break;

		}
	}
	
	// Expansion cutoff and implied param_size ( including (0,0) -> 0 ).
	int M = 1;
	int param_size = 2*M*(M+1);
	
	printf("Path delta set to: %f. \n", path_delta );
	printf("Observation delta set to: %f. \n", obs_delta );
	
	FourierSeries real_v(M);
	real_v.set_mode( 0, 1, 0.5 - 0.5*_Complex_I);
	real_v.set_mode( 1, 0, 0.5 - 0.5*_Complex_I);
	
	Grid real_g( -1, 1, -1, 1, 100, 100 );
	real_g.discretise( &real_v );
	
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
	Vector2d *real_starting_points;
	real_starting_points = (Vector2d*) malloc( parallel_paths*sizeof( Vector2d ) );	
	for( int i = 0; i<parallel_paths; ++i)
		real_starting_points[i] << gsl_rng_uniform(r)*2-1, gsl_rng_uniform(r)*2-1;
	
	Vector2d *observed_starting_points;
	observed_starting_points = (Vector2d*) malloc( parallel_paths*sizeof( Vector2d ) );
	// Generate a path and the observation.
	Vector2d *real_path, *y;
	real_path = ( Vector2d* ) malloc( parallel_paths*path_steps*sizeof( Vector2d ) );
	for( int i=0; i<parallel_paths; ++i )
		em_fast_overdamped_langevin( path_steps, path_delta, C, &real_g, real_path+i*path_steps, r, real_starting_points[i]);

	y = ( Vector2d* ) malloc( parallel_paths*path_steps*sizeof( Vector2d ) );
	for( int i=0; i<parallel_paths; ++i)
	{
		for( int j=0; j<path_steps; ++j )
		{
			noise << gsl_ran_gaussian( r, obs_delta ), gsl_ran_gaussian( r, obs_delta );
			y[i*path_steps+j] = real_path[i*path_steps+j] + noise;
		}
		observed_starting_points[ i ] = y[i*path_steps];
	}
	// Attempt to obtain the value of real_parameter_c with mcmc algo.
	int acceptance_frequency = 0;
	double alpha;
	double acceptance_rate;
	double _Complex average_constant = 0;
	double _Complex *c, *c_star;
	double _Complex *chain = (double _Complex*) malloc(mcmc_trials*param_size*sizeof(double _Complex));
	c = (double _Complex*) malloc(param_size*sizeof(double _Complex));
	c_star = (double _Complex*) malloc(param_size*sizeof(double _Complex));
	double log_u, log_alpha;
	
	// Constants
	for( int j=0; j<param_size; j++ ){
		c[j] = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I;
		//printf("[Debug] c[%i]=%f+%fi.\n",j, creal(c[j]), cimag(c[j]));
		chain[j] = c[j];
	}
	
	double p_param_prob;
	double g_param_prob;
	double p_paths_prob;
		
	FourierSeries v(1);
	v.set_modes( c );
	
	Grid grid( -1, 1, -1, 1, 100, 100 );
	grid.discretise( &v );
	
	Vector2d *x_star = (Vector2d*) malloc( path_steps*parallel_paths*sizeof( Vector2d ) );
	Vector2d *x = (Vector2d*) malloc( path_steps*parallel_paths*sizeof( Vector2d ) );
	
	for( int i=0; i<parallel_paths; ++i )
	{
		printf("[Debug] EM %i \n", i);
		em_fast_overdamped_langevin( path_steps, path_delta, C, &grid, x+i*path_steps, r, observed_starting_points[i]);
	}
	printf("\n ********** \n Starting MCMC. \n ********** \n");
	for( int i=0; i<mcmc_trials+burn; ++i )
	{
		// Propose constants
		for( int j=0; j<param_size; ++j)
			c_star[j] = gsl_ran_gaussian( r, param_delta ) + gsl_ran_gaussian( r, param_delta )*_Complex_I + c[j]; // c* ~ N(c, param_delta).
		
		v.set_modes( c_star );
		grid.discretise( &v );
		for( int i=0; i<parallel_paths; ++i )
			em_fast_overdamped_langevin( path_steps, path_delta, C, &grid, x_star+i*path_steps, r, observed_starting_points[i]); // Generate x*
		
		log_u = log( gsl_rng_uniform(r) );
		p_param_prob = log_p( c_star, param_delta, param_size) - log_p( c, param_delta, param_size);		
		//g_param_prob = log_g( c, c_star, param_delta, param_size) - log_g( c_star, c, param_delta, param_size);
		p_paths_prob = log_p( y, x_star, c_star, path_steps*parallel_paths, obs_delta) - log_p(y, x, c, path_steps*parallel_paths, obs_delta);	
		log_alpha = p_param_prob + g_param_prob + p_paths_prob;
		
		if (log_u<log_alpha)
		{
			printf("[MCMC] n=%i.\n -> Accepted: c_star! \n", i);
			if (i>burn) acceptance_frequency += 1;
			for( int j=0; j<param_size; ++j) c[j] = c_star[j];
			for(int j=0; j<path_steps*parallel_paths; ++j) x[j] = x_star[j];
			v.print_modes();
		}
		// if (i>burn) { // DISABLE BURN IN!
			for( int j=0; j<param_size; ++j)
				chain[ i*param_size + j ] = c[j];
		//}
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
	
	/*
	 * Moving average text file.
	 *
	 */
	double _Complex cumulative_sums[param_size];
	for( int i=0; i<param_size; ++i )
		cumulative_sums[i] = 0.0;
		
	FILE *fp;
	fp = fopen("Testing/MCMC/data_complex.txt", "w+");
	fprintf(fp, "#n");
	
	for(int i=1; i<M+1; i++)
	{
		fprintf(fp,"\tRe(v_\{%i,0\})\tIm(v_\{%i,0\})", i, i);
	}
	
	for(int i=-M; i<M+1; i++)
	for(int j=1; j<M+1; j++)
	{
		fprintf(fp,"\tRe(v_\{%i,%i\})\tIm(v_\{%i,%i\})", i, j, i, j);
	}
	fprintf(fp, "\n");
		
	double _Complex temp_mode;
	int trials_so_far;
	for( int i=0; i<mcmc_trials; i++)
	{
		fprintf(fp, "%i\t", i );
		for( int j=0; j<param_size; ++j)
		{
		cumulative_sums[j] += chain[ i*param_size + j ];
		fprintf(fp, "%f\t%f\t",
					creal(cumulative_sums[j])/(i+1), 
					cimag(cumulative_sums[j])/(i+1) );
		}
		fprintf(fp, "\n");
	}
	fclose(fp);	

	/*
	 * Time series text file.
	 *
	 */
	
	FILE *valf;
	valf = fopen("Testing/MCMC/data_complex_ts.txt", "w+");
	fprintf(valf, "#n");
	
	for(int i=1; i<M+1; i++)
	{
		fprintf(valf,"\tRe(v_\{%i,0\})\tIm(v_\{%i,0\})", i, i);
	}
	
	for(int i=-M; i<M+1; i++)
	for(int j=1; j<M+1; j++)
	{
		fprintf(valf,"\tRe(v_\{%i,%i\})\tIm(v_\{%i,%i\})", i, j, i, j);
	}
	fprintf(valf, "\n");
		
	for( int i=0; i<mcmc_trials; i++)
	{
		fprintf(valf, "%i\t", i );
		for( int j=0; j<param_size; ++j)
		{
			temp_mode = chain[ i*param_size + j ];
			fprintf(valf, "%f\t%f\t", creal(temp_mode), cimag(temp_mode) );
		}
		fprintf(valf, "\n");
	}
	fclose(valf);	
	
	//printf("free(real_path)\n");
	//free(real_path);
	//printf("free(x)\n");
	//free(x);
	//printf("free(x_star)\n");
	//free(x_star);	
	//printf("free(chain)\n");
	//free(chain);
	
}