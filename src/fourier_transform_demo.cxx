#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "fourier_series.h"
#include "grid.h"

int main( int argc, char *argv[] )
{

	int M = 3;
	fftw_complex* in;
	fftw_complex* out;
	in =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*M*M);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*M*M);
	
	printf("in is: \n \n");
	
	for(int i=0; i<M; ++i )
	{
		for(int j=0; j<M; ++j )
		{
			in[i*M+j][0] = 1.0;
			in[i*M+j][1] = 1.0;
		    printf("%f + %fi  \t", in[i*M+j][0], in[i*M+j][1]);
		}
		printf("\n");
	}
	
	fftw_plan plan;
	
	plan = fftw_plan_dft_2d(M, M, in, out, 1, FFTW_ESTIMATE);
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
		
	printf("Result of fft: \n \n");
	for(int i=0; i<M; ++i )
	{
		for(int j=0; j<M; ++j )
			printf("%f + %fi  \t", out[i*M+j][0], out[i*M+j][1]);
		printf("\n");
	}
	
	fftw_free(in);
	fftw_free(out);
		
	double min_x, min_y, max_x, max_y;
	min_x = 0;
	max_x = 1;
	min_y = 0;
	max_y = 1;
	
	const int Nx = 7;
	const int Ny = 7;
	
	FourierSeries f_series(1);
	f_series.set_mode( 1, 1, -0.5*_Complex_I );
	f_series.print_modes();
	
	fftw_complex *y;
	fftw_complex *Y;
	y =  (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	Y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	
	double Lx = max_x-min_x;
	double Ly = max_y-min_y;
	double dx = Lx/Nx;
	double dy = Ly/Ny;
	
	int idx;
	
	printf("f_series.evaluate( 0.25, 0 )=%f. \n", f_series.evaluate( 0.25, 0 ));
	printf("f_series.grad( 0, 0 )=%f. \n", f_series.grad( 0, 0 )(0) );
	
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			y[idx][0] = f_series.evaluate( min_x + i*dx, min_y + j*dy );
			y[idx][1] = 0;
	}
	
	fftw_complex* V;
	fftw_complex* v;
	V = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	v = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny);
	
	fftw_plan plan = fftw_plan_dft_2d(Nx, Ny, y, Y, 1, FFTW_ESTIMATE);
	
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
		
	for( int i=0; i<Nx/2; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			V[idx][0] = -(2*M_PI*i/Lx)*Y[idx][1];
			V[idx][1] =  (2*M_PI*i/Lx)*Y[idx][0];
	}

	for( int i=int(Nx/2); i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			V[idx][0] = -(2*M_PI*(i-Nx)/Lx)*Y[idx][1];
			V[idx][1] =  (2*M_PI*(i-Nx)/Lx)*Y[idx][0];
	}
	
	if (Nx%2==0) {
		printf("[Debug] Nx is even: Nx=%i... \n", Nx );
		for( int j=0; j<Ny; ++j ) {
			idx = (Nx/2)*Ny + j;
			V[idx][0] = 0;
			V[idx][1] = 0;
		}
	}
	
	plan = fftw_plan_dft_2d(Nx, Ny, V, v, -1, FFTW_ESTIMATE);
	if (plan != NULL)
		fftw_execute(plan);
	else
		printf("Aborting: FFTW Plan is NULL.\n");
		
	double padding_coefficient = -Nx*Ny; // Absolutely no idea why this is. Esp. negative.
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			idx = i*Ny + j;
			v[idx][0]/=padding_coefficient;
			v[idx][1]/=padding_coefficient;
			printf("[Debug] v[idx][0]=%f, \n", v[idx][0] );
		}
	
	double x, yy;
	for( int i=0; i<Nx; ++i )
		for( int j=0; j<Ny; ++j ) {
			x = min_x + i*dx;
			yy = min_y + j*dy;
			printf( "Testing d/dx at the point (x,y) = (%f, %f). \n", x, yy );
			printf( "Expected d/dx: %f. Actual d/dx: %f. Ratio: %f \n", 
					f_series.grad(x, yy)(0), v[i*Ny+j][0], v[i*Ny+j][0]/f_series.grad(x, yy)(0) );
		}
	
	fftw_free(Y);
	fftw_free(y);
	fftw_free(V);
	fftw_free(v);

	return 0;
}