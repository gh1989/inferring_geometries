#include <Eigen/Dense>
#include "euler_maruyama.h"
#include <getopt.h>
#include "fourier_series.h"
#include <stdio.h>

using namespace Eigen;

int main(int argc,char *argv[])
{
	bool create_gnu_script = false;
	bool print_data = false;
	char *output_data_file;
	char *gnu_script_file;
	double dt = 0.001;
	double sigma = 0.01;
	int opt;
	int N = 100;
	int seed = 123456789;
	int M = 1;
	int parallel_runs=1;
	gnu_script_file = "Testing/Trajectories/trajectory_plot.gp"; 
	output_data_file = "Testing/Trajectories/data.txt";
	
	// usage / options
	
	/* d - path delta
	 * N - mcmc trials
	 * o - output data file
	 * g - gnuplot script file
	 * G - boolean to create a script file
	 * s - rng seed
	 * M - Fourier mode index cutoff
	 * k - number of parallel running particles
	 * q - diffusion matrix amplitude
	 */
	
	while( ( opt=getopt( argc, argv, ":d:N:o:G:V:s:M:k:q:") )!=EOF ) {
		switch (opt)
		{
			case 'd':
			dt = atof(optarg);
			break;
			
			case 'N':
			N = atoi(optarg);
			break;
			
			case 'o':
			output_data_file = optarg;
			break;
			
			case 'g':
			gnu_script_file = optarg;
			break;
			
			case 'G':
			create_gnu_script = bool(atoi(optarg));
			break;
			
			case 'V':
			print_data = bool(atoi(optarg));
			break;
			
			case 's':
			seed = atoi(optarg);
			break;
			
			case 'M':
			M = atoi(optarg);
			break;
			
			case 'k':
			parallel_runs = atoi(optarg);
			break;
			
			case 'q':
			sigma = atof(optarg);
			break;
		}
	}
	
	// RNG
	const gsl_rng_type *T;
	gsl_ran_discrete_t *g;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	// All trajectories, create N*parallel_runs vectors in an array
	Vector2d *data;
	data = (Vector2d*) malloc(N*parallel_runs*sizeof(Vector2d));
	
	// All random starting points
	Vector2d* starting_points = (Vector2d*) malloc( parallel_runs*sizeof(Vector2d));
	for( int i=0; i<parallel_runs; ++i )
		starting_points[ i ] << gsl_rng_uniform( r )*2-1, gsl_rng_uniform( r )*2-1;
	
	// The potential defined by the following Fourier modes
	FourierSeries v(M);
	v.set_mode( 0, 1, 0.5 - 0.5*_Complex_I);
	v.set_mode( 1, 0, 0.5 - 0.5*_Complex_I);
	
	// The diffusion matrix
	Eigen::Matrix2d C;
	C << sigma, 0, 0, sigma;
	
	for( int i=0; i<parallel_runs; ++i)
		em_overdamped_langevin(N, dt, C, &v, data+i*N, r, starting_points[i]);
	
	if ( print_data )
		for( int j=0; j<N; ++j ){
			for( int i=0; i<parallel_runs; ++i){
				printf("(%f, %f) ", data[i*N+j](0) , data[i*N+j](1) );
				}
			printf("\n");
		}
	
	// Output the data to a text file
	FILE *fp_data;
	fp_data = fopen(output_data_file, "w+");
	fprintf(fp_data, "#t");
	for( int i=0; i<parallel_runs; ++i){
	fprintf(fp_data, "\tX_%i\tY_%i", i, i);
	}
	fprintf(fp_data, "\n");
	
	for( int j=0; j<N; ++j ){
		for( int i=0; i<parallel_runs; ++i){
			fprintf(fp_data, "%i\t%f\t%f", i, data[i*N+j](0), data[i*N+j](1));
		}
		fprintf(fp_data, "\n\n");
	}		
	fclose(fp_data);
	free(data);

	int x, y;
	double _Complex mode;
	if (create_gnu_script)
	{
		FILE *gnu_script;
		gnu_script = fopen(gnu_script_file, "w+");
		fprintf( gnu_script, "set term gif animate delay 3 size 400, 400\n" );
		fprintf( gnu_script, "set output \"trajectory.gif\"\n" );
		fprintf( gnu_script, "f(x,y) = cos(2*pi*x) + sin(2*pi*x) + cos(2*pi*y) + sin(2*pi*y)\n" );
		fprintf( gnu_script, "g(x,y)=0\n" );
		fprintf( gnu_script, "set palette rgbformulae 33,13,10\n" );
		fprintf( gnu_script, "set sample 128\n" );
		fprintf( gnu_script, "set pm3d map\n" );
		fprintf( gnu_script, "set isosample 128\n" );
		fprintf( gnu_script, "set palette defined (-1 \"white\", 0 \"white\", 1 \"red\")\n" );
		fprintf( gnu_script, "do for [n=1:%i] \{\n", N );
		fprintf( gnu_script, "splot [-1:1][-1:1][-1:1] f(x,y) w pm3d,\\\n");
		fprintf( gnu_script, "\"data.txt\" u 2:3:(g($2,$3)) every :::::n w lp t sprintf(\"n=%%i\", n) pt 0,\\\n" );
		
		for( int i=1; i<parallel_runs; ++i )
		{
			x = 2*i+2;
			y = 2*i+3;
			fprintf( gnu_script, "\"data.txt\" u %i:%i:(g($%i,$%i)) every :::::n w lp t sprintf(\"\") pt 0", x, y, x, y );
			if( i<parallel_runs-1 ) fprintf( gnu_script, "," );
			fprintf(gnu_script, "\\\n");
		}
		
		fprintf( gnu_script, "\}");
		fclose( gnu_script );
	}
}
