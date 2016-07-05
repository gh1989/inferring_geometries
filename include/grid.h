#pragma once
#include <Eigen/Dense>
#include "fourier_series.h"
#include <fftw3.h>

using namespace Eigen;

struct Grid
{
	Grid( double, double, double, double, int, int );
	
	double min_x;
	double min_y;
	double max_x;
	double max_y;
	double x_delta;
	double y_delta;
	
	double Lx;
	double Ly;
	
	int size_x;
	int size_y;
	
	Array<double, Dynamic, 1> discrete_domain_x;
	Array<double, Dynamic, 1> discrete_domain_y;
	Array<double, Dynamic, Dynamic> discrete_range;
	Array<double, Dynamic, Dynamic> grad_discrete_range;
	
	Vector2d interpolate_gradient( double x, double y );
	Vector2d interpolate_gradient( Vector2d v );
	double evaluate( double x, double y );
	double evaluate( Vector2d v );
	void populate( FourierSeries *pf_series );
	void calculate_gradient_field();
	
};

Grid::Grid( double min_x_, double min_y_, double max_x_,double  max_y_, int size_x_, int size_y_) : 
		min_x( min_x_ ), min_y( min_y_ ), max_x( max_x_ ), max_y( max_y_ ), size_x(size_x_), size_y(size_y_)
{
	x_delta = (max_x-min_x)/size_x;
	y_delta = (max_y-min_y)/size_y;
	
	Lx = max_x-min_x;
	Ly = max_y-min_y;
	
	discrete_domain_x = Array<double, Dynamic, 1>( size_x );
	discrete_domain_y = Array<double, Dynamic, 1>( size_y );
		
	for( int i=0; i<size_x; ++i )
		discrete_domain_x(i) = min_x + i*x_delta;
	
	for( int j=0; j<size_y; ++j )
		discrete_domain_y(j) = min_y + j*y_delta;
	
	discrete_range = Array<double, Dynamic, Dynamic>( size_x, size_y );
	grad_discrete_range = Array<double, Dynamic, Dynamic>( size_x, 2*size_y );
}


Vector2d Grid::interpolate_gradient( Vector2d v )
{
	return interpolate_gradient( v(0), v(1) );
}

double Grid::evaluate( Vector2d v )
{
	return evaluate( v(0), v(1) );
}

double Grid::evaluate( double x, double y )
{
	return 1.0;
}

Vector2d Grid::interpolate_gradient( double x, double y )
{
	Vector2d ret;
	Vector2d x_diffs;
	Vector2d y_diffs;
	Matrix2d grid_vals;
	
	int lower_x_idx = x / x_delta;
	int upper_x_idx = lower_x_idx + 1;
	
	int lower_y_idx = y / y_delta;
	int upper_y_idx = lower_y_idx + 1;
	
	double x1, x2, y1, y2;
	x1 = discrete_domain_x(lower_x_idx);
	x2 = discrete_domain_x(upper_x_idx);
	y1 = discrete_domain_y(lower_y_idx);
	y2 = discrete_domain_y(upper_y_idx);
	
	double coefficient = 1/(x2-x1)/(y2-y1);
	x_diffs << x2-x, x-x1;
	y_diffs << y2-y, y-y1;
	
	double v11, v12, v21, v22;
	v11 = grad_discrete_range(lower_x_idx, 2*lower_y_idx);
	v21 = grad_discrete_range(lower_x_idx, 2*upper_y_idx);
	v12 = grad_discrete_range(upper_x_idx, 2*lower_y_idx);
	v22 = grad_discrete_range(upper_x_idx, 2*upper_y_idx);
	
	grid_vals << v11, v12, 
				 v21, v22;
	
	ret(0) = coefficient*x_diffs.transpose()*grid_vals*y_diffs;
	
	v11 = grad_discrete_range(lower_x_idx, 2*lower_y_idx+1);
	v21 = grad_discrete_range(lower_x_idx, 2*upper_y_idx+1);
	v12 = grad_discrete_range(upper_x_idx, 2*lower_y_idx+1);
	v22 = grad_discrete_range(upper_x_idx, 2*upper_y_idx+1);
	
	grid_vals << v11, v12, 
				 v21, v22;
	
	ret(1) = coefficient*x_diffs.transpose()*grid_vals*y_diffs;
	
	return ret;
}

void Grid::calculate_gradient_field()
{
	/* Take fft of the discrete range, multiply by ???
	 * then take an ifft and we are left with a discrete
	 * vector field for the gradient of the original 
	 * scalar field.
	 */
	 
	int idx;
	
	fftw_complex* in;
	fftw_complex* out;
	
	in =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
	
	// Is this the best way to get the populated grid?
	for( int i=0; i<size_x; ++i )
	for( int j=0; j<size_y; ++j )
	{
		in[i*size_y+j][0] = discrete_range(i,j);
		in[i*size_y+j][1] = 0.0;
	}

	fftw_plan plan;
	
	plan = fftw_plan_dft_2d( size_x, size_y, in, out, 1, FFTW_ESTIMATE);
	if( plan != NULL )
		fftw_execute(plan);
		
	for( int i=0; i<size_x; ++i )
		for( int j=0; j<size_y; ++j )
		printf("[Debug] out[i*size_y + j ] = %f+%fi. \n", out[i*size_y + j][0],out[i*size_y + j][1] );
		
	fftw_complex *fdx;
	fftw_complex *fdy;
	fftw_complex *dy;
	fftw_complex *dx;

	fdx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
	fdy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
	dx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
	dy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size_x*size_y);
		
	// What has to be multiplied here?
	for( int i=0; i<size_x/2; ++i )
		for( int j=0; j<size_y/2; ++j )
		{
			idx = i*size_y+j;
			fdx[idx][0] = -(2*M_PI*i) / Lx * out[idx][1];
			fdx[idx][1] =  (2*M_PI*i) / Lx * out[idx][0];
			fdy[idx][0] = -(2*M_PI*j) / Ly * out[idx][1];
			fdy[idx][1] =  (2*M_PI*j) / Ly * out[idx][0];
		}
	
	for( int i=size_x/2; i<size_x; ++i )
		for( int j=size_y/2; j<size_y; ++j )
		{
			idx = i*size_y+j;
			fdx[idx][0] = -(2*M_PI*(i-size_x/2)) / Lx * out[idx][1];
			fdx[idx][1] =  (2*M_PI*(i-size_x/2)) / Lx * out[idx][0];
			fdy[idx][0] = -(2*M_PI*(j-size_y/2)) / Ly * out[idx][1];
			fdy[idx][1] =  (2*M_PI*(j-size_y/2)) / Ly * out[idx][0];
		}
	
	if (size_x % 2)
	for( int j=0; j<size_y; ++j ) 
	{
		idx = (size_x/2)*size_y+size_y;
		fdx[idx][0] = 0;
		fdx[idx][1] = 0;
	}
	
	if (size_y % 2)
	for( int i=0; i<size_x; ++i ) 
	{
		idx = i*size_y+size_y/2;
		fdy[idx][0] = 0;
		fdy[idx][1] = 0;
	}
	
	for( int i=0; i<size_x; ++i )
		for( int j=0; j<size_y; ++j )
		{
			printf( "[Debug] fdx[i*size_y + j ] = %f+%fi. \n", fdx[i*size_y+j][0],fdx[i*size_y+j][1] );
		}
	
	plan = fftw_plan_dft_2d(size_x, size_y, fdx, dx, -1, FFTW_ESTIMATE);
	if( plan != NULL )
		fftw_execute(plan);
		
	plan = fftw_plan_dft_2d(size_x, size_y, fdy, dy, -1, FFTW_ESTIMATE);
	if( plan != NULL )
		fftw_execute(plan);
		
	// Is this the best way to repopulate the gradient grid?
	for( int i=0; i<size_x; ++i )
	for( int j=0; j<size_y; ++j )
	{
		grad_discrete_range(i,2*j) = dx[i*size_y+j][0]/(size_x*size_y); 	//Should be real.
		grad_discrete_range(i,2*j+1) = dy[i*size_y+j][0]/(size_x*size_y); 	
	}
	
}

void Grid::populate( FourierSeries *pf_series ) 
{
	Vector2d tmp;

	for( int i=0; i<size_x; ++i )
		for( int j=0; j<size_y; ++j )
		{
			tmp << discrete_domain_x(i), discrete_domain_y(j);
			discrete_range(i,j) =  pf_series->evaluate( tmp );
		}
	
	calculate_gradient_field();
}