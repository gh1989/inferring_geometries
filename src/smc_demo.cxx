#include "fourier_series.h"
#include "smc.h"

void print_tensor( Tensor<double, 4>, int, int, int );
void print_tensor( Tensor<double,2>, int, int);
void print_tensor( Tensor<double, 3>, int, int );

int main(int argc, char *argv[]) 
{
    const int K = 2;  
    const int T = 100;
    const int N = 8;
    const double dt = 0.1;
    const double observation_noise_variance = 0.1;
    const double trajectory_diffusion_sigma = 0.1;
    const int rng_seed = 12345;
  

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
	
	Tensor<double, 2> w(T, N);
    Tensor<double, 4> x(K, T, N, 2);
    Tensor<double, 1> phat(T);
	Tensor<double, 3> y(K, T, 2);
	
    FourierSeries V(1);
    V.set_mode(1, 1, 0.5-0.5*_Complex_I);

    // The observations.
    printf("Here are the observations.\n");
    generate_observations( r, y, K, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    print_tensor(y, K, T);

    // SMC algorithm.
    sequential_monte_carlo( r, x, w, y, phat, V, K, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
	print_tensor(w,T,N);
    print_tensor(x,K,T,N);
    
    return 0;
}

/*
 * Printing helper functions.
 */
 
void print_tensor( Tensor<double,4> x, int K, int T, int N )
{

    for( int k=0; k<K; ++k )
    {
    printf("Independent path: %i... \n", k );
    for( int i=0; i<T; ++i ){
        for( int j=0; j<N; ++j )
            printf("(%.2f,%.2f)\t", x(k,i,j,0), x(k,i,j,1) );
    printf("\n");
    }
    }
}

void print_tensor( Tensor<double,2> w, int T, int N )
{
    for( int i=0; i<T; ++i ){
        for( int j=0; j<N; ++j ) printf("%.6f\t", w(i,j) );
    printf("\n");
    }
}

void print_tensor( Tensor<double,3> y, int K, int T )
{
    for( int k=0; k<K; ++k )
    {
    printf("Independent path: %i... \n ", k);
    for( int i=0; i<T; ++i )
        printf("(%.2f,%.2f)\n", y(k,i,0), y(k,i,1) );
    }
}
