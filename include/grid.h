#pragma once
#include <assert.h> 
#include <Eigen/Dense>
#include "fourier_series.h"
#include <fftw3.h>

using namespace Eigen;

struct Grid
{
	Grid( double, double, double, double, int, int );
	~Grid();

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
	
	int Nx;
	int Ny;
	
	Array<double, Dynamic, 1> discrete_domain_x;
	Array<double, Dynamic, 1> discrete_domain_y;
	Array<double, Dynamic, Dynamic> discrete_range;
	Array<double, Dynamic, Dynamic> grad_discrete_range;
	
	Vector2d interpolate_gradient( double x, double y );
	Vector2d interpolate_gradient( Vector2d v );
	double evaluate( double x, double y );
	double evaluate( Vector2d v );
	void discretise( FourierSeries *pf_series );

	fftw_complex *y;
	fftw_complex *Y;
	fftw_complex *vx;
	fftw_complex *vy;
	fftw_complex *Vx;
	fftw_complex *Vy;
	
};

Grid::Grid( double min_x_, double max_x_, double min_y_ ,double  max_y_, int size_x_, int size_y_) : 
		min_x( min_x_ ), min_y( min_y_ ), max_x( max_x_ ), max_y( max_y_ ), size_x(size_x_), size_y(size_y_)
{
	Nx = size_x;
	Ny = size_y;
	
	x_delta = (max_x-min_x)/size_x;
	y_delta = (max_y-min_y)/size_y;
		
	Lx = max_x-min_x;
	Ly = max_y-min_y;
	
	discrete_domain_x = Array<double, Dynamic, 1>( Nx );
	discrete_domain_y = Array<double, Dynamic, 1>( Ny );
		
	for( int i=0; i<Nx; ++i )
		discrete_domain_x(i) = min_x + i*x_delta;
	
	for( int j=0; j<Ny; ++j )
		discrete_domain_y(j) = min_y + j*y_delta;
	
	discrete_range = Array<double, Dynamic, Dynamic>( Nx, Ny );
	grad_discrete_range = Array<double, Dynamic, Dynamic>( Nx, 2*Ny );	

	y =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	Y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);	
	Vx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	Vy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	vx = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	vy = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
}

Grid::~Grid()
{
	fftw_free(Y);
	fftw_free(y);
	fftw_free(Vx);
	fftw_free(Vy);
	fftw_free(vx);
	fftw_free(vy);
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
	//printf("[Debug] Got here 1...");
	Vector2d ret;
	Vector2d x_diffs;
	Vector2d y_diffs;
	Matrix2d grid_vals;

   	/*
	 * Find the respective index
	 */
	assert( ( x <= max_x ) && ( x >= min_x ) );
	assert( ( y <= max_y ) && ( y >= min_y ) );

	int lower_x_idx = static_cast<int>( (x - min_x ) / x_delta );
	int upper_x_idx = lower_x_idx + 1;
	
	int lower_y_idx = static_cast<int>( (y - min_y) / y_delta );
	int upper_y_idx = lower_y_idx + 1;
	
	double x1, x2, y1, y2;
	
	/*
	if ((lower_x_idx < 0) || (upper_x_idx < 0) || (lower_y_idx < 0) || (upper_y_idx < 0))
		printf("[Problem at: x=%f, y=%f ] lower_x_idx=%i, upper_x_idx=%i, lower_y_idx=%i, upper_y_idx=%i. \n", x, y,
				lower_x_idx, upper_x_idx, lower_y_idx, upper_y_idx );
	
	if ((lower_x_idx >= Nx) || (upper_x_idx >= Nx) || (lower_y_idx >= Ny) || (upper_y_idx >= Ny))
		printf("[Problem at: x=%f, y=%f ] lower_x_idx=%i, upper_x_idx=%i, lower_y_idx=%i, upper_y_idx=%i. \n", x, y,
				lower_x_idx, upper_x_idx, lower_y_idx, upper_y_idx );
	*/

	// These situations all happen on the boundary or close to it.
	if (upper_y_idx > Ny-1)	{
		upper_y_idx--;
		lower_y_idx--;
	}
	
	if (upper_x_idx > Nx-1)	{
		upper_x_idx--;
		lower_x_idx--;
	}
	
	if (lower_x_idx < 0) {
		upper_x_idx++;
		lower_x_idx++;
	}
	
	if (lower_y_idx <0) {
		upper_y_idx++;
		lower_y_idx++;
	}
	
	x1 = discrete_domain_x(lower_x_idx);
	x2 = discrete_domain_x(upper_x_idx);
	y1 = discrete_domain_y(lower_y_idx);
	y2 = discrete_domain_y(upper_y_idx);
	
	double coefficient = 1/(x2-x1)/(y2-y1);
	x_diffs << x2-x, x-x1;
	y_diffs << y2-y, y-y1;
	
	//printf("2...");
	double v11, v12, v21, v22;
  	v11 = grad_discrete_range(lower_x_idx, 2*lower_y_idx);
	v21 = grad_discrete_range(lower_x_idx, 2*upper_y_idx);
	v12 = grad_discrete_range(upper_x_idx, 2*lower_y_idx);
	v22 = grad_discrete_range(upper_x_idx, 2*upper_y_idx);
	
	grid_vals << v11, v12, 
				 v21, v22;
	
	ret(0) = coefficient*x_diffs.transpose()*grid_vals*y_diffs;
	//printf("3...");
	v11 = grad_discrete_range(lower_x_idx, 2*lower_y_idx+1);
	v21 = grad_discrete_range(lower_x_idx, 2*upper_y_idx+1);
	v12 = grad_discrete_range(upper_x_idx, 2*lower_y_idx+1);
	v22 = grad_discrete_range(upper_x_idx, 2*upper_y_idx+1);
	
	grid_vals << v11, v12, 
				 v21, v22;
	//printf("4...");
	ret(1) = coefficient*x_diffs.transpose()*grid_vals*y_diffs;
	//printf("5... \n");
	return ret;
}

void Grid::discretise( FourierSeries *pf_series ) 
{
	Vector2d tmp;
	//printf("[discretise] Got here! \n");
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j )
		{
			tmp << discrete_domain_x(i), discrete_domain_y(j);
			discrete_range(i,j) =  pf_series->evaluate( tmp );
		}
	
	/* Take fft of the discrete range, multiply by 2*pi*k_{1,2}
	 * then take an ifft and we are left with a discrete
	 * vector field for the gradient of the original 
	 * scalar field.
	 */
	 
	int idx;
	//printf("[discretise] Got here! \n");

	
	double Lx = max_x-min_x;
	double Ly = max_y-min_y;
	double dx = Lx/Nx;
	double dy = Ly/Ny;
	//printf("[discretise] Got here! \n");
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			y[idx][0] = discrete_range(i,j);
			y[idx][1] = 0;
	}
	
	//printf("[discretise] Got here! \n");

	fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, y, Y, 1, FFTW_ESTIMATE);
	
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
	
	/*
	 *  For ifft for d/dx. 
	 */
	 
	//printf("[discretise] Got here! \n");
	for( int i=0; i<Nx/2; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			Vx[idx][0] = -(2*M_PI*i/Lx)*Y[idx][1];
			Vx[idx][1] =  (2*M_PI*i/Lx)*Y[idx][0];
	}

	for( int i=int(Nx/2); i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			Vx[idx][0] = -(2*M_PI*(i-Nx)/Lx)*Y[idx][1];
			Vx[idx][1] =  (2*M_PI*(i-Nx)/Lx)*Y[idx][0];
	}
	
	if (Nx%2==0) {
		//printf("[Debug] Nx is even: Nx=%i... \n", Nx );
		for( int j=0; j<Ny; ++j ) {
			idx = (Nx/2)*Ny + j;
			Vx[idx][0] = 0;
			Vx[idx][1] = 0;
		}
	}
	
	/*
	 *  For ifft for d/dy. 
	 */
	 //printf("[discretise] Got here! \n");
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny/2; ++j ) {
			idx = i*Ny + j;
			Vy[idx][0] = -(2*M_PI*j/Lx)*Y[idx][1];
			Vy[idx][1] =  (2*M_PI*j/Lx)*Y[idx][0];
	}

	for( int i=0; i<Nx; ++i )
		for( int j=int(Ny/2); j<Ny; ++j ) {
			idx = i*Ny + j;
			Vy[idx][0] = -(2*M_PI*(j-Ny)/Lx)*Y[idx][1];
			Vy[idx][1] =  (2*M_PI*(j-Ny)/Lx)*Y[idx][0];
	}
	
	if (Ny%2==0) {
		//printf("[Debug] Ny is even: Nx=%i... \n", Nx );
		for( int i=0; i<Nx; ++i ) {
			idx = i*Ny + Ny/2;
			Vy[idx][0] = 0;
			Vy[idx][1] = 0;
		}
	}
	
	
	plan = fftw_plan_dft_2d(Nx, Ny, Vx, vx, -1, FFTW_ESTIMATE);
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
		
	double padding_coefficient = -Nx*Ny; // Absolutely no idea why this is. Esp. negative.
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			vx[idx][0]/=padding_coefficient;
			vx[idx][1]/=padding_coefficient;
			//printf("[Debug] vx[idx][0]=%f, \n", vx[idx][0] );
		}
	
	plan = fftw_plan_dft_2d(Nx, Ny, Vy, vy, -1, FFTW_ESTIMATE);
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
	
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			vy[idx][0]/=padding_coefficient;
			vy[idx][1]/=padding_coefficient;
			//printf("[Debug] vy[idx][0]=%f, \n", vy[idx][0] );
		}
		

		
	// Is this the best way to repopulate the gradient grid?
	for( int i=0; i<Nx; ++i )
	for( int j=0; j<Ny; ++j )
	{
		grad_discrete_range(i,2*j) = vx[i*Ny+j][0]; 	//Should be real.
		grad_discrete_range(i,2*j+1) = vy[i*Ny+j][0]; 	
	}
}
