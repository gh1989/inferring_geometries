#include "smc.h"

void generate_observations( gsl_rng *r, 
                            Tensor<double, 3> &y, 
                            int K, int T, FourierSeries &V, double dt, 
                            double observation_noise_variance, double trajectory_diffusion_sigma )
{
   /*
    * Generate a path and some observations. X_t follows Langevin SDE.
    * Then the observation is the path at time t plus Gaussian noise with 
    * variance observation_noise_variance.
    */

    Vector2d x;
    Vector2d diffusion_noise;
    
    for( int k=0; k<K; ++k )
    {
        x << gsl_rng_uniform(r)*2-1,gsl_rng_uniform(r)*2-1;
        for( int t=1; t<T; ++t )
        {
            y(k, t, 0) = x(0) + gsl_ran_gaussian(r, observation_noise_variance );
            y(k, t, 1) = x(1) + gsl_ran_gaussian(r, observation_noise_variance );

            diffusion_noise << gsl_ran_gaussian( r, dt ),
                               gsl_ran_gaussian( r, dt );

            x += -V.grad(x)*dt + trajectory_diffusion_sigma*diffusion_noise;
        }
    } 
}


void generate_particle_samples( gsl_rng *r, 
                                Tensor<double, 4> &x, 
                                Tensor<double, 3> &y,
                                int K, int N, int t, FourierSeries &V, double dt, 
                                double observation_noise_variance, double trajectory_diffusion_sigma )
{
    if (t == 0)
    {
        for( int k=0; k<K; ++k )
        for( int i=0; i<N; ++i ){
            x(k,0,i,0) = y(k,0,0) + gsl_ran_gaussian(r, observation_noise_variance);
            x(k,0,i,1) = y(k,0,1) + gsl_ran_gaussian(r, observation_noise_variance);
        }    
        return;    
    }

    Vector2d dx(0,0);
    Vector2d current;
    Vector2d diffusion_noise(0,0);    
    
    // Sample by forward simulation of the previous resampled x_{t-1}.
    // x_t ~ q(.|x_{t-1}, y_t)
    for( int k=0; k<K; ++k )
    for( int i=0; i<N; ++i ){
        current << x(k,t-1,i,0), x(k,t-1,i,0);
        diffusion_noise <<  gsl_ran_gaussian( r, dt ),
                            gsl_ran_gaussian( r, dt );
     
         
        dx = -V.grad(current)*dt + trajectory_diffusion_sigma*diffusion_noise;
        x(k,t,i,0) = x(k,t-1,i,0) + dx(0);
        x(k,t,i,1) = x(k,t-1,i,1) + dx(1);
    }
}


double assign_weights( Tensor<double,4> &x, 
                       Tensor<double,2> &w,
                       Tensor<double,3> &y,
                       int K, int N, int t,
                       double observation_noise_variance )
{
   /*
    * Calculating and storing the unnormalised weights. The sampler used for
    * resampling does not require normalisation for the weights, and those 
    * used for the estimate for the marginal density for the MH acceptance
    * ratio uses unnormalised weights anyway.
    *
    * Each element is a 2d vector, the weight is a bivariate normal with mean
    * x^i_t and variance the identity multiplied by observation_noise_variance.
    */
    
    double exponent_constant = 0.5 / observation_noise_variance;
    double coefficient_constant = 0.5 / (M_PI*observation_noise_variance);
    double total = 0;
    double x_square;
    double y_square;
    
    if( t==0 )
    {
        for( int i=0; i<N; ++i )
        {
            w(t,i) = 1.0/N;
        }
    return 1.0;
    }

    for( int i=0; i<N; ++i )
    {
        x_square = 0;
        y_square = 0;

        for( int k=0; k<K; ++k )
        {
            x_square += (y(k, t, 0)-x(k, t, i, 0))*(y(k, t, 0)-x(k, t, i, 0));
            y_square += (y(k, t, 1)-x(k, t, i, 1))*(y(k, t, 1)-x(k, t, i, 1));
        }

        w(t,i) = coefficient_constant*exp(-exponent_constant*(x_square+y_square));
    }
    
    for( int i=0; i<N; ++i )
        total += w(t,i);
    for( int i=0; i<N; ++i )
        w(t,i) /= total;
    return total;
}

double estimate_marginal_likelihood( int t, int N, Tensor<double,2> &w, Tensor<double, 1> &phat )
{
   /* 
    * Estimate of the marginal likelihood. 
    * \hat{p}_c(y_0) = 1/N \sum_{ i=1 }^N w_0^{*i}.
    * \hat{p}_c( { \bf y_t } ) = 1/N \hat{p}_c( {\bf y_{t-1} } ) * \sum_{i=1}^N w_t^{*i}
    * ( The unnormalised weights are summed. )
    */
    double total = 0;
    for( int i=0; i<N; ++i )
        total += w(t, i);
    
    if( t==0 )
        return (1.0/N) * total;

    return (1.0/N) * phat( t-1 ) * total;
}

double sequential_monte_carlo(  gsl_rng *r,
                                Tensor<double,4> &x,
                                Tensor<double,2> &w,
                                Tensor<double,3> &y,
                                Tensor<double, 1> &phat,
                                FourierSeries &V,
                                int K, int N, int T, double dt, 
                                double observation_noise_variance, double trajectory_diffusion_sigma )
{
   /*
    * Start by generating samples and assign weights at time t=0, calculate the
    * estimated marginal likelihood for time 0. In the main loop resample with
    * replacement N times given the weights (which do not need to be normalised
    * according to GSL documentation on gsl_ran_discrete function usage.)
    */

    //printf("Trying sequential monte carlo with potential modes: \n");
    //V.print_modes();
    
    generate_particle_samples(r, x, y, K, N, 0, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    double norm_total;    
    norm_total = assign_weights( x, w, y, K, N, 0, observation_noise_variance );
    phat(0) = norm_total*estimate_marginal_likelihood(0, N, w, phat);
    
    Tensor<double, 4> resampled( K, T, N, 2 );
    gsl_ran_discrete_t *g;

    unsigned int resample_particle_index;
    double log_estimated_marginal;
    double p[N];
    
    for( int t=1; t<T; t++ )
    {
        // Now to resample the particles with replacement.
        for( int i=0; i<N; ++i )
            p[i] = w(t,i);
        g = gsl_ran_discrete_preproc(N, p);
    
        // Resampling: work out which to replace.
        for(int i=0; i<N; ++i ){
            // Sample {0, ..., N-1} with weights given by p.
            resample_particle_index = gsl_ran_discrete(r, g);
            
            for( int k=0; k<K; ++k )
                for( int j=0; j<t; ++j ){
                    resampled(k, j,i,0) = x(k, j, resample_particle_index, 0 );
                    resampled(k, j,i,1) = x(k, j, resample_particle_index, 1 );
            }        
        }
        
        // Resampling: replace.
        for(int i=0; i<N; ++i ){
            for( int j=0; j<t; ++j ){
                for( int k=0; k<K; ++k ){
                    x(k, j, i, 0) = resampled(k, j,i,0);
                    x(k, j, i, 1) = resampled(k, j,i,1);
                }
            }
        }
                  
        generate_particle_samples(r, x, y, K, N, t, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
        norm_total = assign_weights( x, w, y, K, N, t, observation_noise_variance );
        phat(t) = norm_total*estimate_marginal_likelihood(t, N, w, phat);
    }
    
    log_estimated_marginal = log( phat(T-1) );
    return log_estimated_marginal;
}
