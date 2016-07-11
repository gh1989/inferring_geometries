#include "fourier_series.h"
#include "mcmc.h"

int main(int argc, char *argv[]) 
{
	// Settings
	int mcmc_trials = 100;
	int parallel_paths = 12;
	double param_noise = 0.05;
	double path_noise = 0.1;
	double path_delta = 0.01;
	int path_length = 100;
	int rng_seed = 1234;
	double sigma = 0.01;
	
	const int M = 1;
	FourierSeries v(M);
	v.set_mode(1,1,0.5);

	Matrix2d C;
	C << sigma, 0, 0, sigma;

	Vector2d noise;

	// Trajectory containers
	// The true paths
	Vector2d *real_x;
	real_x = (Vector2d*) malloc( sizeof(Vector2d) * path_length * parallel_paths * mcmc_trials );

	// The observed paths
	Vector2d *y;
	y = (Vector2d*) malloc( sizeof(Vector2d) * path_length * parallel_paths * mcmc_trials );
 
	// The starting positions of each particle (sample path)
	Vector2d *real_starting_points;
	real_starting_points = (Vector2d*) malloc( sizeof(Vector2d) * parallel_paths );

     	// Setup the random number generator
	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set( r, rng_seed );

	for( int i=0; i<parallel_paths; ++i )
	{
		real_starting_points[i] << gsl_rng_uniform(r)*2 - 1, gsl_rng_uniform(r)*2 - 1;
	}

	int idx;
	// Generate observations
	for( int k=0; k<parallel_paths; ++k ){
	    em_overdamped_langevin( path_length, path_delta, C, &v, real_x+k*path_length, r, real_starting_points[i]); 
	    for( int j=0; j<path_length; ++j ){
		noise << gsl_ran_gaussian(r, obs_noise), gsl_ran_gaussian(r, obs_noise);
		idx = k*path_length + j;
		y[idx] = real_x[idx] + noise;
		}
	}

	double _Complex parameter;
	parameter = gsl_ran_gaussian(r, param_noise) + _Complex_I* gsl_ran_gaussian(r, param_noise);
	double sigma_param;
	sigma_param = gsl_ran_lognormal(r, 0, param_noise );	

	Vector2d *x;
	x = (Vector2d*) malloc( sizeof(Vector2d)*path_length*parallel_paths*mcmc_trials );

	for( int k=0; k<parallel_paths; ++k){
		x[k*path_length] = y[k*path_length] + gsl_rng_gaussian(r, sigma_param);
	 }

	double *weights;
	weights = (double*)malloc( sizeof(double)*parallel_paths*path_length );

	for( int k=0; k<parallel_paths; ++k )
	{
		weights[
	}

	return 0;
}
