#include "fourier_series.h"
#include "smc.h"

void print_tensor( Tensor<double, 3>, int, int );
void print_tensor( Tensor<double,2>, int, int);
void print_tensor( Tensor<double, 2>, int );

int main(int argc, char *argv[]) 
{
    const int N = 8;
    const int T = 100;
    const double dt = 0.1;
    const double observation_noise_variance = 0.1;
    const double trajectory_diffusion_sigma = 0.1;
    const int rng_seed = 12345;
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
	
	Tensor<double, 2> w(T, N);
    Tensor<double, 3> x(T, N, 2);
    Tensor<double, 1> phat(T);
	Tensor<double, 2> y(T, 2);
	
    FourierSeries V(1);
    V.set_mode(1, 1, 0.5);

    // The observations.
    printf("Here are the observations.\n");
    generate_observations( r, y, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    print_tensor(y, T);

    // SMC algorithm.
    sequential_monte_carlo( r, x, w, y, phat, V, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
	print_tensor(w,T,N);
    print_tensor(x,T,N);
    
    return 0;
}

/*
 * Printing helper functions.
 */
 
void print_tensor( Tensor<double,3> x, int r, int c )
{
    for( int i=0; i<r; ++i ){
        for( int j=0; j<c; ++j )
            printf("(%.2f,%.2f)\t", x(i,j,0), x(i,j,1) );
    printf("\n");
    }
}

void print_tensor( Tensor<double,2> x, int r, int c )
{
    for( int i=0; i<r; ++i ){
        for( int j=0; j<c; ++j ) printf("%.6f\t", x(i,j) );
    printf("\n");
    }
}

void print_tensor( Tensor<double,2> x, int l )
{
    for( int i=0; i<l; ++i )
        printf("(%.2f,%.2f)\n", x(i,0), x(i,1) );
}
