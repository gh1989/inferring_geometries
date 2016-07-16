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
								int N, int t )
{
	
	if (t == 0)
	{
		for( int i=0; i<N; ++i )
		{
			// Initial density.
			x(0,i,0) = gsl_rng_uniform(r);
			x(0,i,1) = gsl_rng_uniform(r);
		}	
		return;	
	}

}


void assign_weights( Tensor<double,3> &x, 
					 Tensor<double,2> &w,
					 Tensor<double,2> &y,
					 int N, int t,
					 double observation_noise_variance )
{
	if (t == 0)
	{
		for( int i=0; i<N; ++i )
		{
			// Initial density.
			w(0,i) = calculate_weight( y, x, N, i, 0, observation_noise_variance );
		}	
		return;	
	}

}

double calculate_weight( Tensor<double,2> &y, Tensor<double,3> &x, 
						 int N, int i, int t,
						 double observation_noise_variance )
{
	double normalisation=0;
	double weight=0;
	double exponent_constant = 	(0.5 / observation_noise_variance);
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

	printf("Weight calculated(t=%i, i=%i): %f.\n", t, i, weight );

	return weight;
}

