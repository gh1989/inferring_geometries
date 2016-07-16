#include "fourier_series.h"
#include "smc.h"

void print_tensor( Tensor<double, 3>, int, int );
void print_tensor( Tensor<double,2>, int, int);
void print_tensor( Tensor<double, 2>, int );

int main(int argc, char *argv[]) 
{
	// Number of particles.
	const int N = 4;

	// Time horizon.
	const int T = 8;

	// Prepare random number generator (rng).
	const int rng_seed = 12345;
	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set( r, rng_seed );

	// The path timestep.
	const double dt = 0.1;
	const double observation_noise_variance = 0.1;
	const double trajectory_diffusion_sigma = 0.1;
	

	Tensor<double, 2> w(T, N);
	Tensor<double, 3> x(T, N, 2);

	// The projection of the potential onto the finite-d subspace
	FourierSeries V(1);
	V.set_mode(1, 1, 0.5);

	// The observations.
	printf("Here are the observations.\n");
	Tensor<double,2> y(T, 2);
	generate_observations( r, y, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
	print_tensor(y, T);

	// Initialising.
	generate_particle_samples(r, x, N, 0, V, dt, trajectory_diffusion_sigma );
	// The particles.
	printf("Here are the particles.\n");
	// A quick test...
	print_tensor( x, T, N);

	// The weights.
	assign_weights( x, w, y, N, 0, observation_noise_variance );
	printf("Here are the weights.\n");
	// A quick test...
	print_tensor( w, T, N);

	//Calculation of the marginal likelihood estimate
	double phat_total = 0;
	Tensor<double,1> phat(T);
	for( int i=0; i<N; ++i )
		phat_total += w( 0, i );
	phat_total /= N;
	phat(0) = phat_total;
	printf("We found that the marginal likelihood estimate for t=0: %f. \n", phat(0) ); 

	Tensor<double, 3> resampled( T, N, 2 );
	gsl_ran_discrete_t *g;

	unsigned int resample_particle_index;
	double p[N];
	
	for( int t=1; t<T; t++ )
	{
		
		// Now to resample the particles with replacement, weighted with their
		// respective weights as calculated earlier.
		for( int i=0; i<N; ++i )
			p[i] = w(t,i);

		// Resample with replacement moving forward.
		g = gsl_ran_discrete_preproc(N, p);
	
		for(int i=0; i<N; ++i )
		{
			resample_particle_index = gsl_ran_discrete(r, g);
			for( int j=0; j<t; ++j )
			{
				resampled(j,i,0) = x( j, resample_particle_index, 0 );
				resampled(j,i,1) = x( j, resample_particle_index, 1 );
			}		
		}

		// Now resample in place...
		for(int i=0; i<N; ++i )
		{
			for( int j=0; j<t; ++j )
			{
				x(j, i, 0) = resampled(j,i,0);
				x(j, i, 1) = resampled(j,i,1);
			}
		}
		
		printf("Resampling the particles with replacement. (t=%i) \n", t);	
		// Check what the particles look like now we have resampled with replacement.
		print_tensor(x,T,N);

		printf("Moving forward. (t=%i) \n", t);
		// Move all x forward
		generate_particle_samples(r, x, N, t, V, dt, trajectory_diffusion_sigma );
		print_tensor(x,T,N);

		// Recalculate the weights.
		assign_weights( x, w, y, N, t, observation_noise_variance );
		printf("Recalculating weights. (t=%i) \n", t);
		print_tensor(w,T,N);

		//Calculation of the marginal likelihood estimate
		phat_total = 0;
		
		for( int i=0; i<N; ++i )
		{		
			phat_total += w( t, i );
		}
		phat_total /= N; 
		phat(t) = phat(t-1)*phat_total;
	}
	return 0;
}

void print_tensor( Tensor<double,3> x, int r, int c )
{
	for( int i=0; i<r; ++i )
	{
		for( int j=0; j<c; ++j )
		{
		printf("(%.2f,%.2f)\t", x(i,j,0), x(i,j,1) );
		}
	printf("\n");
	}
}

void print_tensor( Tensor<double,2> x, int r, int c )
{
	for( int i=0; i<r; ++i )
	{
		for( int j=0; j<c; ++j )
		{
		printf("%.6f\t", x(i,j) );
		}
	printf("\n");
	}
}

void print_tensor( Tensor<double,2> x, int l )
{
	for( int i=0; i<l; ++i )
	{
	printf("(%.2f,%.2f)\n", x(i,0), x(i,1) );
	}
}
