#include "fourier_series.h"
const int M = 1;

int main(){

	FourierSeries V(M);
	V.set_mode( 1, 0, 0.5-0.5*_Complex_I);
	V.set_mode( 0, 1, 0.5-0.5*_Complex_I);
	
	double _Complex temp_mode;
	for( int i=-M; i<M+1; ++i)
	for( int j= -M; j<M+1; ++j )
	{
		temp_mode = V.get_mode( i, j );
		printf("[Debug] V.get_mode( %i, %i ) = %f+%f i \n", i, j, creal( temp_mode ), cimag( temp_mode ) );
	}
	
	Vector2d v(0, 0);
	Vector2d dv = V.grad( v );
	
	printf("[Debug] gradient_field (grad V)(%f,%f) = (%f,%f). \n", v(0), v(1), dv(0), dv(1) );

	
	v(0)=1;
	v(1)=1;
	dv = V.grad( v );
	
	printf("[Debug] gradient_field (grad V)(%f,%f) = (%f,%f). \n",v(0), v(1), dv(0), dv(1) );
	
	v(0)=0.25;
	v(1)=0.25;
	dv = V.grad( v );
	
	printf("[Debug] gradient_field (grad V)(%f,%f) = (%f,%f). \n",v(0), v(1), dv(0), dv(1) );
	
	double _Complex parameters[2*M*(M+1)];
	
	int ii;
	int jj;
	int offset;
	int idx;
	
	/* Attempt to generate the following: (0,0) is always 0. 
	 *
	 *          2 3 4    2 3 4
	 *			. . 1 -> 1 0 1
	 *			. . .    4 3 2
	 */	
	
	// First half row
	for( int i=1; i<M+1; ++i )
		parameters[(i-1)] = i+i*_Complex_I;
	
	// Remaining parts of the upper plane
	for( int j=1; j<M+1; ++j )
	for( int i=-M; i<M+1; ++i )
	{
		int kk = 2*M+(j-1)*(2*M+1)+i+1;
		printf("[Debug] %i\n", kk);
		parameters[ 2*M+(j-1)*(2*M+1)+i ] = kk+_Complex_I*kk;
	}
	
	V.set_modes( parameters );
	V.print_modes();
	
	return 0;
}