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
	const int rng_seed = 1234;
	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set( r, rng_seed );

	// The path timestep.
	const double dt = 0.01;
	const double observation_noise_variance = 0.1;
	const double trajectory_diffusion_sigma = 0.01;
	

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
	generate_particle_samples( r, x, N, 0 );

	// The particles.

	printf("Here are the particles.\n");
	// A quick test...
	print_tensor( x, T, N);

	assign_weights( x, w, y, N, 0, observation_noise_variance );
	// The weights.
	printf("Here are the weights.\n");
	// A quick test...
	print_tensor( w, T, N);

	//Calculation of the marginal likelihood estimate
	double phat_0 = 0;
	Tensor<double,1> phat(T);
	for( int i=0; i<N; ++i )
		phat_0 += w( 0, i );
	phat_0 /= N;
	phat(0) = phat_0;
	printf("We found that the marginal likelihood estimate for t=0: %f. \n", phat(0) ); 

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
