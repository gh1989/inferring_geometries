#include <stdio.h>
#include "mcmc.h"
#include <Eigen/Dense>
#include <getopt.h>
#include <gsl/gsl_cdf.h>
#include <math.h>

using namespace Eigen;

int main(int argc, char *argv[] )
{
	int opt=0;

	// S: MAP SIZE.
	int map_size = 1024;
	// N: TRIALS.
	int N = 1;
	// T: GENERATIONS
	int T = 128;
	// B: BURN PROBABILITY
	double burn_chance = 0.45;

	Array<int, Dynamic, Dynamic> terrain( map_size, map_size );
	Array<int, Dynamic, Dynamic> current( map_size, map_size );	
	Array<int, Dynamic, Dynamic> distribution( map_size, map_size );
	
	double u;
	int rng_seed = 123;
	int infection_direction_x;
	int infection_direction_y;
	int infect_x;
	int infect_y;
	char file_name[38];
	// Setup the random number generator
	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc( gsl_rng_default );
	gsl_rng_set( r, rng_seed );

	FILE *fp;

	for( int x=0; x<map_size; ++x )
		for( int y=0; y<map_size; ++y )
		{
			terrain(x,y) = 0;
			distribution(x,y)=0;
		}
	// Put a fire in the middle.
	terrain( map_size/2, map_size/2 ) = 1;	

	for( int i=0; i<N; ++i )
	{
		for( int x=0; x<map_size; ++x )
			for( int y=0; y<map_size; ++y )
			{
				if( terrain( x, y ) == 1 )
				current( x, y ) = 1;
				else
				current( x, y ) = 0;
			}
		
		for( int t=0; t<T; ++t )
		{
			for( int x=0; x<map_size; ++x )
				for( int y=0; y<map_size; ++y )
					if( current(x,y) == 1 )
					{
						u = gsl_rng_uniform(r);
						if( u < burn_chance )
						{
							infection_direction_x = gsl_rng_uniform_int(r, 3)-1;
							infection_direction_y = gsl_rng_uniform_int(r, 3)-1;
							infect_x = x+infection_direction_x;
							infect_y = y+infection_direction_y;
							if( infect_x < 0 ) infect_x = 0;
							if( infect_x > map_size - 1 ) infect_x = map_size - 1;
							if( infect_y < 0 ) infect_y = 0;
							if( infect_y > map_size - 1 ) infect_y = map_size - 1;
							printf(" (%i, %i) burns (%i,%i). \n",x,y, infect_x, infect_y);
							current(infect_x, infect_y) = 1;
						}
					}

		sprintf(file_name, "burn_frames/burn_distribution_%04u.txt", t );
		printf("Got here... \n");
		fp = fopen(file_name, "w+");
		for( int x=0; x<map_size; x++ )
		{
		for( int y=0; y<map_size; y++ )
			fprintf(fp,"%i\t", current(x,y));
		fprintf(fp,"\n");	
		}	
		fclose(fp);
		}
	for( int x=0; x<map_size; x++ )
	for( int y=0; y<map_size; y++ )
		distribution(x,y) += current(x,y);
	

	}
	
	return 0;
}

