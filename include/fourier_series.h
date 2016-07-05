#pragma once
#include <complex.h>
#undef I
#include <Eigen/Dense>
#include <math.h>

using namespace Eigen;

// M is the cutoff for |i| <= M.
// A full set -M <= i = M is 2*M+1 values
// Due to the reality constraint v*_{-k} = v_{k}.
const double _Complex two_pi_i = 2*M_PI*_Complex_I;

struct FourierSeries
{
	Vector2d grad( Vector2d );
	double evaluate( Vector2d );
	void set_mode( int,  int, double _Complex );
	void set_modes( double _Complex* );
	double _Complex get_mode( int, int );
	void get_modes( double _Complex* );
	FourierSeries(int M_);
	~FourierSeries();
	double _Complex *modes;
	int total_modes;
	void print_modes();
	int M;
};

FourierSeries::FourierSeries( int M_ ) : M(M_), total_modes( (2*M_+1)*(2*M_+1) )
{
	//printf("[Debug] FourierSeries(%i) \n", M_ );
	modes = (double _Complex*) malloc( total_modes*sizeof( double _Complex ) );
	for( int i = 0; i<total_modes; ++i )
	    modes[i] = 0;
}

FourierSeries::~FourierSeries()
{
	//printf("[Debug] ~FourierSeries() \n" );
	free(modes); // Seg fault... why?
}
		
void FourierSeries::set_mode( int i, int j, double _Complex mode)
{
	int idx;
	int ii = M + i;
	int jj = M + j;
	idx = ii + jj*(2*M+1);
	
	modes[idx] = mode;
	//printf("[Debug] set_mode: modes[%i] = %f + %fi. \n", idx, creal(modes[idx]), cimag(modes[idx]) );
	
	// The reality constraint:
	ii = M - i;
	jj = M - j;
	idx = ii + jj*(2*M+1);
	modes[idx] = conj(mode);
	//printf("[Debug] set_mode: modes[%i] = %f + %fi. \n", idx, creal(modes[idx]), cimag(modes[idx]) );
	
	//print_modes();
}

void FourierSeries::get_modes( double _Complex* c )
{
	int ii, jj, idx;
	for(int i=-M; i<M+1; ++i)
	for(int j= 0; j<M+1; ++j)
	{
		ii = i + M;
		idx = ii + j*(M+1);
		c[ idx ] = get_mode(i, j);
	}
}

void FourierSeries::set_modes( double _Complex* c)
{
	/* Attempt to generate the following: (0,0) is always 0. 
	 *
	 *   7  8  9 10 11		7  8  9 10 11
	 *   2  3  4  5  6		2  3  4  5  6
	 *   .  .  .  0  1  => 	1  0  0  0  1
	 *   .  .  .  .  .		6  5  4  3  2
	 *   .  .  .  .  . 		11 10 9  8  7 
	 */

	int ii, jj;
	int idx_setting_from;
	int offset;
	double _Complex tmp;
	
	// First half row  ( 0, 1 above )
	for( int i=1; i<M+1; ++i )
	{
		tmp = c[(i-1)];
		//printf("[Debug] set_modes: (i,j)=(%i,%i). c[%i] = %f + %fi. \n", i, 0, idx, creal(tmp), cimag(tmp) );
		set_mode( i, 0, tmp );
	}
	
	// Remaining parts of the upper plane (2, ..., 11 above )
	for( int j=1; j<M+1; ++j )
	for( int i=-M; i<M+1; ++i )
	{
		idx_setting_from = 2*M+(j-1)*(2*M+1)+i;
		tmp = c[idx_setting_from];
		set_mode( i, j, tmp );
		//printf("[Debug] set_modes: (i,j)=(%i,%i). c[%i] = %f + i %f. \n", i, j, idx_setting_from, creal(tmp), cimag(tmp) );
	}
	
}

double _Complex FourierSeries::get_mode( int i, int j )
{
	if (i==0 && j==0) 
		return 0;
	int ii = M + i;
	int jj = M + j;
	int idx = ii + jj*(2*M+1);
	double _Complex the_mode = modes[idx];
	// printf("[Debug] get_mode: modes[%i] = %f + %fi. \n", idx, creal(modes[idx]), cimag(modes[idx]) );
	return the_mode;
}

Vector2d FourierSeries::grad( Vector2d v )
{	
	Vector2d ret;
	
	double x = v(0);
	double y = v(1);
	
	double _Complex x_component;
	double _Complex y_component;	
	double _Complex mode;
	double _Complex tmp;
	
	for( int i=-M; i<M+1; ++i)
	for( int j= -M; j<M+1; ++j )
	{
		// Get the mode for (i, j) and compute the terms for the Fourier series.
		// One is for the positive, one is for the negative values of j which is
		// given by the reality constraint.
		mode = get_mode( i, j );
		
		// printf("[Debug] Vector2d FourierSeries<%i> grad(v): (%i,%i). \n", M, i, j );
		
		tmp = mode * cexp(  two_pi_i*(i*x+j*y) );	
		x_component += i*two_pi_i*tmp;
		y_component += j*two_pi_i*tmp;
		
		
		/*printf("[Debug] x_component (%f,%f). y_component (%f,%f). \n", 
				creal( x_component ), cimag( x_component), 
				creal( y_component),  cimag( y_component) );	
		
		*/
	}
	
	ret(0) = creal( x_component );
	ret(1) = creal( y_component );
	
	return ret;
}

double FourierSeries::evaluate( Vector2d v )
{
	double _Complex ret = 0;
	double x = v(0);
	double y = v(1);
	
	for( int i=-M; i<M+1; ++i )
		for( int j=-M; j<M+1; ++j )
				ret += cexp(  two_pi_i*(i*x+j*y) );
				
	return creal( ret );
}

/*
 * Helper functions.
 */

void FourierSeries::print_modes()
{
	printf("Potential modes: \n");
	double _Complex md;
	for( int j=M; j>-M-1; --j )
	{
		for( int i=-M; i<M+1; ++i )
		{
			md = get_mode(i,j);
			printf(" v_\{%i, %i\} = %f + %fi. \t", i, j, creal(md), cimag(md));
		}
		printf("\n");
	}
}