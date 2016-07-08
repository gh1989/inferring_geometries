#include "euler_maruyama.h"



void em_ornstein_uhlenbeck( int n, double dt, double c, Vector2d *path_generated, gsl_rng *r )
{
	// Starting point
	Vector2d x;
	x << gsl_rng_uniform( r )*2-1, gsl_rng_uniform(r)*2-1;
	
	double root_dt = sqrt( 2*dt );
	Vector2d eps;
	path_generated[0] = x;
	for( int i=1; i<n; ++i )
	{
		eps << gsl_ran_gaussian(r,1), gsl_ran_gaussian(r,1);
		x += -c*x*dt + root_dt*eps;
		path_generated[ i ] = x;
	}
}

int em_overdamped_langevin(int n, double dt, Matrix2d C, FourierSeries *V, Vector2d *data, gsl_rng *r, Vector2d start_from)
{
	Vector2d xi;
	Vector2d x;
	data[0] = start_from;
	double rootdt = sqrt( 2*dt );
	x << start_from(0), start_from(1);
	
	for( int i = 1; i<n; ++i )
	{
		xi << gsl_ran_gaussian(r,1), gsl_ran_gaussian(r,1); 
		x += -V->grad( x )*dt + C*xi*rootdt;
		data[i] = x;
	}
	
	return 0;
}

int em_fast_overdamped_langevin( int n, double dt, Matrix2d C, Grid *V, Vector2d *data, gsl_rng *r, Vector2d start_from)
{
	Vector2d xi;
	Vector2d x;
	
	double min_x = V->min_x;
	double max_x = V->max_x;
	double min_y = V->min_y;
	double max_y = V->max_y;
	
	double temp_x;
	double temp_y;
	
	data[0] = start_from;
	double rootdt = sqrt( 2*dt );
	x << start_from(0), start_from(1);
	
	for( int i = 1; i<n; ++i )
	{
		xi << gsl_ran_gaussian(r,1), gsl_ran_gaussian(r,1); 

		// Assuming that x is never so far out of bounds that the offset
		// would not be enough.
		if( x(0) < min_x ) 
			x(0) += (max_x-min_x);
		if( x(0) > max_x ) 
			x(0) -= (max_x-min_x);
		if( x(1) < min_y ) 
		    x(1) += (max_y-min_y);
		if( x(1) > max_y ) 
			x(1) -= (max_y-min_y);

		x += -V->interpolate_gradient( x )*dt + C*xi*rootdt;

		// Assuming that x is never so far out of bounds that the offset
		// would not be enough.
		if( x(0) < min_x ) 
			x(0) += (max_x-min_x);
		if( x(0) > max_x ) 
			x(0) -= (max_x-min_x);
		if( x(1) < min_y ) 
		    x(1) += (max_y-min_y);
		if( x(1) > max_y ) 
			x(1) -= (max_y-min_y);
				
		data[i] = x;
	}
	
	return 0;
}
