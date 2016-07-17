#include "mcmc.h"

// single constant
double log_p( double c, double param_delta)
{
    /* p( c ) = N( 0, delta )
     * p( c ) = 1/2pi e^( -0.5/delta*c^2 )
     * log p( c ) =* -0.5/delta*c^2 upto a constant.
     */
    return -0.5/param_delta*c*c;
}

double log_p( double _Complex c, double param_delta )
{
    return -0.5/param_delta*cabs(c)*cabs(c);
}

double log_p_sigma( double sigma, double param_delta )
{
    return -0.5/param_delta*( log(sigma) )*( log(sigma) );
}

double log_g_sigma( double sigma_star, double sigma, double param_delta)
{
    return -0.5/param_delta*( log(sigma_star)-log(sigma) )*( log(sigma)-log(sigma) );
}

// sequence of complex constants
double log_p( double _Complex *c, double param_delta, int length )
{
    /* p( c ) = N( 0, delta )
     * p( c ) = 1/2pi e^( -0.5/delta*c^2 )
     * log p( c ) =* -0.5/delta*c^2 upto a constant.
     */
     
    double sum = 0.0;
     
    for( int i=0; i<length; ++i )
        sum += -(0.5/param_delta)*cabs(c[i])*cabs(c[i]);
     
    return sum;
}

double log_g( double c_star, double c, double param_delta )
{
    /* g( c* | c ) = N( c, delta )
     * this means that g( c* | c ) = 1/2pi e^( -0.5/delta*(c-c*)^2 )
     * log g( c* | c ) =* -0.5*param_delta(c-c*)^2 upto a constant.
     */
    return -0.5/param_delta*(c_star-c)*(c_star-c);
}

double log_g( double _Complex *c_star, double _Complex *c, double param_delta, int length )
{
    /* g( c* | c ) = N( c, delta )
     * this means that g( c* | c ) = 1/2pi e^( -0.5/delta*(c-c*)^2 )
     * log g( c* | c ) =* -0.5*param_delta(c-c*)^2 upto a constant.
     */
     
     double sum = 0.0;
     for (int i=0; i<length; ++i )
        sum += cabs( c_star[i] - c[i] )*cabs( c_star[i] - c[i] );
        
    return -0.5/param_delta*sum;
}

double log_p( Vector2d *y, Vector2d *x, int path_steps, double obs_delta, int real_per_observed)
{
    /* p( y | x, c ) for the Ornstein Uhlenbeck process
     * y_i = x_i + eps
     * so that y_i ~ N( x_i, delta ) since eps ~ N( 0, delta ).
     * p( y_i | x_i, c ) = 1/sqrt(2*pi*delta)*e^( -0.5/delta*(y_i-x_i)^2 )
     * They are conditionally independent so that:
     * p( y_0, ... , y_n | x_0, ..., x_n, c ) = p( y_0 | x_0, c )...p( y_n | x_n, c )
     * log p( y | x, c ) = p( y_0 | x_0, c ) + ... + p(y_n | x_n, c )
     *                     =* -0.5/obs_delta*(x_0-y_0)^2 - ... -0.5/obs_delta*(x_n-y_n)^2
     * upto a constant, which vanishes in the acceptance probability calculation.
     */

    double log_p_sum = 0.0;
    for( int i=0; i<path_steps; ++i )
        log_p_sum -= (0.5/obs_delta)*(y[i] - x[i*real_per_observed]).transpose()*(y[i] - x[i*real_per_observed]);
    return log_p_sum;
}