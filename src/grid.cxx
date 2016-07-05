#include "grid.h"


int main(){
	
	Array<double, Dynamic, 1> discrete_domain(128);
	
	double max_x, max_y, min_x, min_y;
	
	max_x = 1.0;
	max_y = 1.0;
	min_x = 0.0;
	min_y = 0.0;
	
	const int nx = 100;
	const int ny = 100;
	
	Grid discrete_function_approximation(min_x, min_y, max_x, max_y, nx, ny );
	
	/*
	for( int i=0; i<nx; ++i )
		printf("[Grid] x[%i]=%f. \n", i, discrete_function_approximation.discrete_domain_x(i));
		
	for( int i=0; i<ny; ++i )
		printf("[Grid] y[%i]=%f. \n", i, discrete_function_approximation.discrete_domain_y(i));
	*/
	
	FourierSeries f_series(1);
	f_series.set_mode( 1, 1, 0.5 - 0.5*_Complex_I );
	discrete_function_approximation.populate( &f_series );
	
	/*
	for( int i=0; i<nx; ++i )
	for( int j=0; j<ny; ++j )
	{
		printf("[Grid] grid(%i,%i)=%f. ", i, 2*j, discrete_function_approximation.discrete_range(i,2*j));
		printf("[Grid] grid(%i,%i)=%f. \n", i, 2*j+1, discrete_function_approximation.discrete_range(i,2*j+1));
	}
	*/
	Vector2d expected, output;
	double x,y;
	
	x=0.0;
	y=0.0;
	output = discrete_function_approximation.interpolate_gradient( x, y );
	expected = f_series.grad( Vector2d(x,y) );
	printf("[Grid] approximation \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, output(0), output(1)); 
	printf("[Grid] expected \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, expected(0), expected(1)); 
	
	x=0.125;
	y=0.25;
	output = discrete_function_approximation.interpolate_gradient( x, y );
	expected = f_series.grad( Vector2d(x,y) );
	printf("[Grid] approximation \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, output(0), output(1)); 
	printf("[Grid] expected \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, expected(0), expected(1)); 
		
	x=0.06125;
	y=0.0125;
	output = discrete_function_approximation.interpolate_gradient( x, y );
	expected = f_series.grad( Vector2d(x,y) );
	printf("[Grid] approximation \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, output(0), output(1)); 
	printf("[Grid] expected \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, expected(0), expected(1)); 
	
	x=0.030625;
	y=0.06125;
	output = discrete_function_approximation.interpolate_gradient( x, y );
	expected = f_series.grad( Vector2d(x,y) );
	printf("[Grid] approximation \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, output(0), output(1)); 
	printf("[Grid] expected \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, expected(0), expected(1)); 
	
	return 0;
}