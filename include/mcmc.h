#pragma once
#include <Eigen/Dense>
#include "euler_maruyama.h"
#include "fourier_series.h"
#include <gsl/gsl_cdf.h>
#include <math.h>
#include <stdio.h>

// single constant
double log_p( double c, double param_delta)
{
	/* p( c ) = N( 0, delta )
	 * p( c ) = 1/2pi e^( -0.5/delta*c^2 )
	 * log p( c ) =* -0.5/delta*c^2 upto a constant.
	 */
	return -0.5/param_delta*c*c;
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

double log_p( Vector2d *y, Vector2d *x, double c, int path_steps, double obs_delta)
{
	/* p( y | x, c ) for the Ornstein Uhlenbeck process
	 * y_i = x_i + eps
	 * so that y_i ~ N( x_i, delta ) since eps ~ N( 0, delta ).
	 * p( y_i | x_i, c ) = 1/sqrt(2*pi*delta)*e^( -0.5/delta*(y_i-x_i)^2 )
	 * They are conditionally independent so that:
	 * p( y_0, ... , y_n | x_0, ..., x_n, c ) = p( y_0 | x_0, c )...p( y_n | x_n, c )
	 * log p( y | x, c ) = p( y_0 | x_0, c ) + ... + p(y_n | x_n, c )
	 *					 =* -0.5/obs_delta*(x_0-y_0)^2 - ... -0.5/obs_delta*(x_n-y_n)^2
	 * upto a constant, which vanishes in the acceptance probability calculation.
	 */

	double log_p_sum = 0.0;
	for( int i=0; i<path_steps; ++i )
		log_p_sum -= (0.5/obs_delta)*(y[i] - x[i]).transpose()*(y[i] - x[i]);
	return log_p_sum;
}

double log_p(  Vector2d *y, Vector2d *x, double _Complex *c, int path_steps, double obs_delta)
{
	double hack_c = 1.0; // c is not explicitly used.
	return log_p(  y, x, hack_c, path_steps, obs_delta);
}