#include "fourier_series.h"
#include "smc.h"
#include <getopt.h>
#include "util.h"

int main(int argc, char *argv[]) 
{
    int K = 1;  
    int T = 16;
    int N = 8;
    int M = 8;
    double ds = 0.001;
    double observation_noise_variance = 0.1;
    double trajectory_diffusion_sigma = 0.1;
    int rng_seed = 123;
  
    int opt;

    while( ( opt = getopt( argc, argv, ":T:R:K:N:L:d:o:t:c:M:" ) ) != EOF ) 
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
            case 'M':
            M = atoi(optarg);
            break;
            case 'd':
            ds = atof(optarg);
            break;
            case 'o':
            observation_noise_variance = atof(optarg);
            break;  
            case 't':
            trajectory_diffusion_sigma = atof(optarg);
            break;   
        }
    }
    
    double dt = M*ds;
    int augmented_total = (T-1)*M+T;
    
    gsl_rng *r;
    gsl_rng_env_setup();
    r = gsl_rng_alloc( gsl_rng_default );
    gsl_rng_set( r, rng_seed );
	
	Tensor<double, 2> w(augmented_total, N);
    Tensor<double, 4> x(K, augmented_total, N, 2);
    Tensor<double, 1> phat(augmented_total);
	Tensor<double, 3> y(K, T, 2);
    Tensor<double, 4> resampled( K, augmented_total, N, 2 );
    
    FourierSeries V(1);
    V.set_mode(1, 1, 0.5-0.5*_Complex_I);
    
    // The observations.
    printf("Here are the observations.\n");
    generate_observations( r, y, K, T, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    print_tensor(y, K, T);
    
    for( int t=0; t<T-1; t++ )
    {
    generate_particle_samples( r, 
                               x, y,
                               K, N, t, M,
                               V, 
                               dt, ds,
                               observation_noise_variance, 
                               trajectory_diffusion_sigma );
    }
    print_tensor(x, K, augmented_total, N );
    
    return 0;
}
