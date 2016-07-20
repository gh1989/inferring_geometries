#include "fourier_series.h"
#include "smc.h"
#include <getopt.h>
#include "util.h"

int main(int argc, char *argv[]) 
{
    int K = 1;  
    int T = 16;
    int N = 8;
    double dt = 0.001;
    double observation_noise_variance = 0.1;
    double trajectory_diffusion_sigma = 0.1;
    int rng_seed = 123;
  
    int opt;

    while( ( opt = getopt( argc, argv, ":T:R:K:N:L:d:o:t:c:" ) ) != EOF ) 
    {
    switch (opt)
        {
            case 'T':
            T = atoi(optarg);
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
            case 'd':
            dt = atof(optarg);
            break;
            case 'o':
            observation_noise_variance = atof(optarg);
            break;  
            case 't':
            trajectory_diffusion_sigma = atof(optarg);
            break;   
        }
    }

    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
	
	Tensor<double, 2> w(T, N);
    Tensor<double, 4> x(K, T, N, 2);
    Tensor<double, 1> phat(T);
	Tensor<double, 3> y(K, T, 2);
    Tensor<double, 4> resampled( K, T, N, 2 );
    
    FourierSeries V(1);
    V.set_mode(1, 1, 0.5-0.5*_Complex_I);

    // The observations.
    printf("Here are the observations.\n");
    generate_observations( r, y, K, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    print_tensor(y, K, T);

    // SMC algorithm.
    sequential_monte_carlo( r, x, w, y, phat, resampled, V, K, N, T, dt, observation_noise_variance, trajectory_diffusion_sigma );
	print_tensor(w,T,N);
    print_tensor(x,K,T,N);
    
    return 0;
}
