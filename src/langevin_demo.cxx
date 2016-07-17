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
    int path_steps = 100;
    int seed = 123456789;
    int M = 1;
    int parallel_paths=1;
    gnu_script_file = "Testing/Trajectories/trajectory_plot.gp"; 
    output_data_file = "Testing/Trajectories/data.txt";
    
    // usage / options
    
    /* d - path delta
     * N - path steps
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
            path_steps = atoi(optarg);
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
            parallel_paths = atoi(optarg);
            break;
            
            case 'q':
            sigma = atof(optarg);
            break;
        }
    }
    
    printf("path delta: %f \n", dt);
    printf("path steps: %i \n", path_steps);
    printf("output data file: %s\n", output_data_file);
    printf("gnuplot script file: %s \n", gnu_script_file);
    printf("create a script file: %i \n", create_gnu_script);
    printf("rng seed: %i \n", seed);
    printf("Fourier mode index cutoff: %i\n", M);
    printf("number of parallel running particles: %i\n", parallel_paths);
    printf("diffusion matrix amplitude: %f\n",sigma);
        
    
    // RNG
    const gsl_rng_type *T;
    gsl_ran_discrete_t *g;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // Trajectory create path_steps*parallel_paths vectors in an array
    Vector2d *data;
    data = ( Vector2d* ) malloc( parallel_paths*path_steps*sizeof( Vector2d ) );
    
    // All random starting points
    Vector2d* starting_points;
    starting_points = (Vector2d*) malloc( parallel_paths*sizeof( Vector2d ) );
    
    for( int i=0; i<parallel_paths; ++i )
        starting_points[i] << 2*gsl_rng_uniform( r )-1, 2*gsl_rng_uniform( r )-1;
        
    // The potential defined by the following Fourier modes
    FourierSeries v(M);
    v.set_mode( 0, 1, 0.5 - 0.5*_Complex_I);
    v.set_mode( 1, 0, 0.5 - 0.5*_Complex_I);
    
    // The diffusion matrix
    Eigen::Matrix2d C;
    C << sigma, 0, 0, sigma;
            
    for( int i=0; i<parallel_paths; ++i)
    {
        em_overdamped_langevin(path_steps, dt, C, &v, data+i*path_steps, r, starting_points[i]);
    }
    
    if ( print_data )
        for( int j=0; j<path_steps; ++j ){
            for( int i=0; i<parallel_paths; ++i){
                printf("(%f, %f) ", data[i*path_steps+j](0) , data[i*path_steps+j](1) );
                }
            printf("\n");
        }

    // Output the data to a text file
    FILE *fp_data;
    fp_data = fopen(output_data_file, "w+");
    fprintf(fp_data, "#t");
    for( int i=0; i<parallel_paths; ++i){
    fprintf(fp_data, "\tX_%i\tY_%i", i, i);
    }
    fprintf(fp_data, "\n");
    
    for( int j=0; j<path_steps; ++j ){
        for( int i=0; i<parallel_paths; ++i){
            fprintf(fp_data, "%i\t%f\t%f", i, data[i*path_steps+j](0), data[i*path_steps+j](1));
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
        fprintf( gnu_script, "do for [n=1:%i] {\n", path_steps );
        fprintf( gnu_script, "splot [-1:1][-1:1][-1:1] f(x,y) w pm3d,\\\n");
        fprintf( gnu_script, "\"data.txt\" u 2:3:(g($2,$3)) every :::::n w lp t sprintf(\"n=%%i\", n) pt 0,\\\n" );
        
        for( int i=1; i<parallel_paths; ++i )
        {
            x = 2*i+2;
            y = 2*i+3;
            fprintf( gnu_script, "\"data.txt\" u %i:%i:(g($%i,$%i)) every :::::n w lp t sprintf(\"\") pt 0", x, y, x, y );
            if( i<parallel_paths-1 ) fprintf( gnu_script, "," );
            fprintf(gnu_script, "\\\n");
        }
        
        fprintf( gnu_script, "}");
        fclose( gnu_script );
    }
}
