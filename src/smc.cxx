#include "smc.h"

void generate_observations( gsl_rng *r, 
                            Tensor<double, 2> &y, 
                            int T, FourierSeries &V, double dt, 
                            double observation_noise_variance, double trajectory_diffusion_sigma )
{

    Vector2d x(0,0);
    Vector2d diffusion_noise(0,0);

    for( int t=1; t<T; ++t )
    {
        y(t,0) = x(0) + gsl_ran_gaussian(r, observation_noise_variance );
        y(t,1) = x(1) + gsl_ran_gaussian(r, observation_noise_variance );

        diffusion_noise << gsl_ran_gaussian( r, 2*dt ),
                           gsl_ran_gaussian( r, 2*dt );

        x += -V.grad(x)*dt + trajectory_diffusion_sigma*diffusion_noise;
    } 
}


void generate_particle_samples( gsl_rng *r, 
                                Tensor<double, 3> &x, 
                                Tensor<double, 2> &y,
                                int N, int t, FourierSeries &V, double dt, 
                                double observation_noise_variance, double trajectory_diffusion_sigma )
{
    
    if (t == 0)
    {
        for( int i=0; i<N; ++i ){
            x(0,i,0) = y(0,0) + gsl_ran_gaussian(r, observation_noise_variance);
            x(0,i,1) = y(0,1) + gsl_ran_gaussian(r, observation_noise_variance);
        }    
        return;    
    }

    Vector2d dx(0,0);
    Vector2d current;
    Vector2d diffusion_noise(0,0);    
    
    // Sample by forward simulation of the previous resampled x_{t-1}.
    // x_t ~ q(.|x_{t-1}, y_t)
    for( int i=0; i<N; ++i ){
        current << x(t-1,i,0), x(t-1,i,0);
        diffusion_noise << gsl_ran_gaussian( r, 2*dt ),
                         gsl_ran_gaussian( r, 2*dt );
     
         
        dx = -V.grad(current)*dt + trajectory_diffusion_sigma*diffusion_noise;
        x(t,i,0) = x(t-1,i,0) + dx(0);
        x(t,i,1) = x(t-1,i,1) + dx(1);
    }
}


void assign_weights( Tensor<double,3> &x, 
                     Tensor<double,2> &w,
                     Tensor<double,2> &y,
                     int N, int t,
                     double observation_noise_variance )
{
    for( int i=0; i<N; ++i )
        w(t,i) = calculate_weight( y, x, N, i, t, observation_noise_variance );  
}

double calculate_weight( Tensor<double,2> &y, Tensor<double,3> &x, 
                         int N, int i, int t,
                         double observation_noise_variance )
{
    double normalisation=0;
    double weight=0;
    double exponent_constant =     (0.5 / observation_noise_variance);
    double x_square=0;
    double y_square=0;

    for( int j=0; j<N; ++j )
    {
        x_square = (y(t,0)-x(t,j,0))*(y(t,0)-x(t,j,0));
        y_square = (y(t,1)-x(t,j,1))*(y(t,1)-x(t,j,1));
        normalisation += exp(-exponent_constant*(x_square+y_square));
    }

    x_square = (y(t,0)-x(t,i,0))*(y(t,0)-x(t,i,0));
    y_square = (y(t,1)-x(t,i,1))*(y(t,1)-x(t,i,1));
    weight = exp(-exponent_constant*(x_square+y_square))/normalisation;
    return weight;
}

void sequential_monte_carlo( gsl_rng *r,
                             Tensor<double,3> &x,
                             Tensor<double,2> &w,
                             Tensor<double,2> &y,
                             FourierSeries &V,
                             int N, int T, double dt, 
                             double observation_noise_variance, double trajectory_diffusion_sigma )
{
    generate_particle_samples(r, x, y, N, 0, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
    assign_weights( x, w, y, N, 0, observation_noise_variance );

    Tensor<double, 3> resampled( T, N, 2 );
    gsl_ran_discrete_t *g;

    unsigned int resample_particle_index;
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
            
            for( int j=0; j<t; ++j ){
                resampled(j,i,0) = x( j, resample_particle_index, 0 );
                resampled(j,i,1) = x( j, resample_particle_index, 1 );
            }        
        }
        
        // Resampling: replace.
        for(int i=0; i<N; ++i ){
            for( int j=0; j<t; ++j ){
                x(j, i, 0) = resampled(j,i,0);
                x(j, i, 1) = resampled(j,i,1);
            }
        }
                  
        generate_particle_samples(r, x, y, N, t, V, dt, observation_noise_variance, trajectory_diffusion_sigma );
        assign_weights( x, w, y, N, t, observation_noise_variance );
     		
    }
}