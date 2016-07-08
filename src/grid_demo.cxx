#include "grid.h"

void bounds_test( double x, double y, double min_x, double max_x, double min_y, double max_y, double x_delta, double y_delta);;
void index_test( double v, double min_v, double max_v, double v_delta, int size_v );
void interpolation_gradient_test( double x, double y, Grid *g, FourierSeries *f );

int main(){
	double max_x, max_y, min_x, min_y;
	
	max_x = 1.0;
	max_y = 1.0;
	min_x = -1.0;
	min_y = -1.0;
	
	const int nx = 100;
	const int ny = 100;
	
	double x_delta = ( max_x - min_x ) / nx;
	double y_delta = ( max_y - min_y ) / ny;
	
	printf("[Debug] Expect: 0.5, 0.5. \n");
	bounds_test( 0.5, 0.5, min_x, max_x, min_y, max_y, x_delta, y_delta);
	printf("[Debug] Expect: -0.25, 0.5. \n");
	bounds_test( -0.25, 0.5, min_x, max_x, min_y, max_y, x_delta, y_delta);;
	printf("[Debug] Expect: 0.25, -0.5. \n");
	bounds_test( 0.25, 1.5, min_x, max_x, min_y, max_y, x_delta, y_delta);
	printf("[Debug] Expect: 0.75, 0.5. \n");
	bounds_test( -1.25, 0.5, min_x, max_x, min_y, max_y, x_delta, y_delta);
	
	
	printf("[Debug] Expect: lower_idx:0 , upper_idx:1. \n");
	index_test(-0.995, min_x, max_x, x_delta, nx);
	printf("[Debug] Expect: lower_idx:1 , upper_idx:2. \n");
	index_test(-0.979, min_x, max_x, x_delta, nx);
	printf("[Debug] Expect: lower_idx:98, upper_idx:99. \n");
	index_test(0.99, min_x, max_x, x_delta, nx);
	printf("[Debug] Expect: lower_idx:98 , upper_idx:99. \n");
	index_test(1.0, min_x, max_x, x_delta, nx);
	
	printf("[Debug] Got here! \n");
	Grid discrete_function_approximation(min_x, max_x, min_y, max_y, nx, ny );
	printf("[Debug] Got here! \n");
	FourierSeries f_series(1);
	printf("[Debug] Got here! \n");
	f_series.set_mode( 1, 1, -0.5*_Complex_I );
	printf("[Debug] Got here! \n");
	discrete_function_approximation.discretise( &f_series );
	printf("[Debug] Got here! \n");
	interpolation_gradient_test(0.0, 0.0, &discrete_function_approximation, &f_series);
	interpolation_gradient_test(0.125, 0.25, &discrete_function_approximation, &f_series);
	interpolation_gradient_test(0.06125, 0.0125, &discrete_function_approximation, &f_series);
	interpolation_gradient_test(0.030625, 0.06125, &discrete_function_approximation, &f_series);
	// It is not the responsibility of the interpolation on periodic BCs to ensure that the vector is within bounds
	interpolation_gradient_test(-1.0, 0.06125, &discrete_function_approximation, &f_series);	
	//interpolation_gradient_test(1.0, 4.06125, &discrete_function_approximation, &f_series);
	
	return 0;
}

void bounds_test( double x, double y, double min_x, double max_x, double min_y, double max_y, 
				  double x_delta, double y_delta )
{
	// while would work more robustly, no while loops...
	if( x < min_x ) 
		x += (max_x-min_x);
	
	if( x > max_x ) 
		x -= (max_x-min_x);
	
	// while would work more robustly, no while loops...
	if( y < min_y ) 
		y += (max_y-min_y);
	
	if( y > max_y ) 
		y -= (max_y-min_y);
		
	printf("[Debug] Result: %f, %f. \n", x , y);
	
}

void index_test( double v, double min_v, double max_v, double v_delta, int size_v )
{
	int lower_idx = static_cast<int>( (v - min_v) / v_delta );
	int upper_idx = lower_idx + 1;
	
	while (upper_idx > size_v-1)
	{
		lower_idx--;
		upper_idx--;
	}
	
	while (lower_idx < 0)
	{
		lower_idx++;
		upper_idx++;
	}
	
	printf("[Debug] lower_idx: %i, upper_idx:%i. \n * \n", lower_idx, upper_idx);
}

void interpolation_gradient_test( double x, double y, Grid *g, FourierSeries *f )
{
	Vector2d expected, output;
	output = g->interpolate_gradient( x, y );
	expected = f->grad( Vector2d(x,y) );
	printf("[Debug] approximation \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, output(0), output(1)); 
	printf("[Debug] expected \tgrad( %f, %f ) =\t (%f, %f)\n", x, y, expected(0), expected(1)); 
	
}
