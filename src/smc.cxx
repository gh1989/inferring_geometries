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
        // Each one starts in a random [-1,1]^2 grid.
        x << gsl_rng_uniform(r)*2-1, gsl_rng_uniform(r)*2-1;
        y(k, 0, 0) = x(0) + gsl_ran_gaussian(r, observation_noise_variance );
        y(k, 0, 1) = x(1) + gsl_ran_gaussian(r, observation_noise_variance );
        
        for( int t=1; t<T; ++t )
        {
            // Now transition x forwards according to Langevin.
            diffusion_noise << gsl_ran_gaussian( r, 1 ), gsl_ran_gaussian( r, 1 );
            x += -V.grad(x)*dt + trajectory_diffusion_sigma*sqrt(dt)*diffusion_noise;

            // The vector y_t is a noisy measurement of x_t at time t.
            y(k, t, 0) = x(0) + gsl_ran_gaussian(r, observation_noise_variance );
            y(k, t, 1) = x(1) + gsl_ran_gaussian(r, observation_noise_variance );
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
            x(k,0,i,0) = y(k,0,0); //+ gsl_ran_gaussian(r, observation_noise_variance);
            x(k,0,i,1) = y(k,0,1); //+ gsl_ran_gaussian(r, observation_noise_variance);
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
        current << x(k,t-1,i,0), x(k,t-1,i,1);
        diffusion_noise <<  gsl_ran_gaussian( r, 1 ), gsl_ran_gaussian( r, 1 );
        dx = -V.grad(current)*dt + trajectory_diffusion_sigma*sqrt(dt)*diffusion_noise;
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
    double diff_total; 
    for( int i=0; i<N; ++i )
    {
        diff_total = 0;
        for( int k=0; k<K; ++k )
        {
            diff_total += (y(k, t, 0)-x(k, t, i, 0))*(y(k, t, 0)-x(k, t, i, 0));
            diff_total += (y(k, t, 1)-x(k, t, i, 1))*(y(k, t, 1)-x(k, t, i, 1));
        }
    w(t,i) = -exponent_constant*diff_total;
    }

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
        total += exp( w(t, i) );
    
    if( t==0 )
        return (1.0/N) * total;

    return (1.0/N) * phat( t-1 ) * total;
}

double sequential_monte_carlo( gsl_rng *r,
                               Tensor<double, 4> &x,
                               Tensor<double, 2> &w,
                               Tensor<double, 3> &y,
                               Tensor<double, 1> &phat,
                               Tensor<double, 4> &resampled,
                               FourierSeries &V,
                               int K, int N, int T, double dt, 
                               double observation_noise_variance,
                               double trajectory_diffusion_sigma )
{
   /*
    * Start by generating samples and assign weights at time t=0, calculate the
    * estimated marginal likelihood for time 0. In the main loop resample with
    * replacement N times given the weights (which do not need to be normalised
    * according to GSL documentation on gsl_ran_discrete function usage.)
    */
    gsl_ran_discrete_t *g;
    
    for( int t=0; t<T; t++ )
    {
        generate_particle_samples( r, x, y, K, N, t, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
        assign_weights( x, w, y, K, N, t, observation_noise_variance );
        phat(t) = estimate_marginal_likelihood(t, N, w, phat);
        resample( r, g, x, w, resampled, K, N, t );
    }
    
    return log( phat(T-1) );
}

void resample( gsl_rng *r,
               gsl_ran_discrete_t *g,
               Tensor<double,4> &x,
               Tensor<double,2> &w,
               Tensor<double, 4> &resampled,
               int K, int N, int t )
{
    unsigned int resample_particle_index;
    double p[N];
    
    for( int i=0; i<N; ++i )
        p[i] = exp( w(t, i) );
    g = gsl_ran_discrete_preproc(N, p);  
    
    for( int i=0; i<N; ++i )
    {
        resample_particle_index = gsl_ran_discrete(r, g);
        for( int k=0; k<K; ++k )
        {
            resampled(k, t, i, 0) = x(k, t, resample_particle_index, 0 );
            resampled(k, t, i, 1) = x(k, t, resample_particle_index, 1 );
        }
    }
    
    // Resampling: replace.
    for(int i=0; i<N; ++i )
    for( int k=0; k<K; ++k ){
        x(k, t, i, 0) = resampled(k, t, i, 0);
        x(k, t, i, 1) = resampled(k, t, i, 1);
    }
}