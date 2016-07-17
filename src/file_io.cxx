#include "file_io.h"

void output_sigma_chain( int mcmc_trials, double *sigma_chain )
{
    FILE *fp;
    fp = fopen("Testing/MCMC/data_sigma_ts.txt", "w+");    
    
    for( int i=0; i<mcmc_trials; ++i )
    {
        fprintf(fp, "%f\n", sigma_chain[i]);
    }
    fclose(fp);
    
    fp = fopen("Testing/MCMC/data_sigma_average_ts.txt", "w+");
    double running_average = 0;
    
    for( int i=0; i<mcmc_trials; ++i )
    {
        running_average+=sigma_chain[i];
        fprintf(fp, "%f\n", running_average/(i+1));
    }
    fclose(fp);
}

void output_posterior_density_complex( int M, int mcmc_trials, int param_size, double _Complex *chain, int bins){
   /*
    * Output for the histogram. This function takes the chain, which contains the parameter chains for all 
    * parameters (:param_size: many of them) continguous in memory. So much of the logic is ensuring that this 
    * data is unpacked in the right order. 
    *
    * This function works out the minimum and maximum values (real and complex) for each mode in the chain,
    * splits it into :bins: intervals and counts number of elements in each interval in order to output the
    * results in a histogram or heatmap in gnuplot.
    */
    
    // Create a pointer for the output file.
    FILE *fp;
    fp = fopen("Testing/MCMC/data_posterior_density_complex.txt", "w+");    
        
    int idx;
    
    double _Complex tmp_value;
    double min_real, min_imag, max_real, max_imag;
    double real_delta, imag_delta;
    double tmp_real, tmp_imag;
    double tmp_range_imag, tmp_range_real;
    double tmp_range_imag_next, tmp_range_real_next;
    
    double *range_x = (double*) malloc( sizeof(double)*bins*param_size ); 
    double *range_y = (double*) malloc( sizeof(double)*bins*param_size ); 
    double *weights = (double*) malloc( sizeof(double)*bins*bins*param_size );

    // For every parameter
    for( int p=0; p<param_size; ++p )
    {
        
        // Get the first elements of each parameter chain, and let these be min/max/imag/real.
        min_real = creal(chain[p]);
        min_imag = cimag(chain[p]);
        max_real = creal(chain[p]);
        max_imag = cimag(chain[p]);
    
        // Find the max/min/imag/real
        for( int i=0; i<mcmc_trials; ++i )
        {
            idx = i*param_size+p; // pth parameter on the ith trial.
            tmp_value = chain[idx];
            tmp_real = creal(tmp_value);
            tmp_imag = cimag(tmp_value);
            if( tmp_real < min_real ) min_real = tmp_real;
            if( tmp_real >= max_real ) max_real = tmp_real;
            if( tmp_imag < min_imag ) min_imag = tmp_imag;
            if( tmp_imag >= max_imag ) max_imag = tmp_imag;
        }
    
        // Now split the interval up into the bins
        real_delta = ( max_real - min_real ) / bins;
        imag_delta = ( max_imag - min_imag ) / bins; 
    
        // The real (x) and imaginary (y) dimensions and points
        for( int i=0; i<bins; ++i )
        {
            range_x[p*bins+i] = min_real + (0.5 + i ) * real_delta;
            range_y[p*bins+i] = min_imag + (0.5 + i ) * imag_delta;
        }
    
        // Prepare the weights
        for( int i=0; i<bins; ++i )
        for( int j=0; j<bins; ++j )
            weights[p*bins*bins+i*bins+j] = 0;

        // Now go through the whole chain for a given parameter and count how many times it falls within
        // each of the intervals for the histogram.
        for( int i=0; i<mcmc_trials; ++i )
        {
            idx = i*param_size+p;
            tmp_value = chain[idx];
            tmp_imag = cimag( tmp_value );
            tmp_real = creal( tmp_value );
            
            for( int j=0; j<bins; ++j ){
                tmp_range_real = range_x[p*bins+j];
                tmp_range_real_next = (j==bins-1)?tmp_range_real:range_x[p*bins+j+1];
                
                if( (  tmp_real >= tmp_range_real ) && ( tmp_real <= tmp_range_real_next ) ) {
                    for( int k=0; k<bins; ++k ) {
                    
                        tmp_range_imag = range_y[p*bins+k];
                        tmp_range_imag_next = (k==bins-1)? tmp_range_imag:range_y[p*bins+k+1];
                        if( ( tmp_imag >= tmp_range_imag ) && ( tmp_imag <= tmp_range_imag_next ) ) {
                            weights[p*bins*bins+j*bins+k]+=1.0;
                            break;
                        }
                    }
                    break;
                }
            }            
            
        }

    }

    // Finally write the file.
    for( int i=0; i<bins; ++i )
    for( int j=0; j<bins; ++j )
    {
        for( int p=0; p<param_size; ++p ) 
        {
            tmp_range_imag = range_y[p*bins+j];
            tmp_range_real = range_x[p*bins+i];
            fprintf( fp, "%f\t%f\t%f\t", tmp_range_real, tmp_range_imag,  weights[p*bins*bins+i*bins+j] );
        }
        fprintf(fp, "\n");
    }
    
    free( range_x );
    free( range_y );
    fclose(fp);
}

void output_time_series_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain ){
    /*
     * Time series text file.
     *
     */
    double _Complex temp_mode;
    FILE *fp;
    fp = fopen("Testing/MCMC/data_complex_ts.txt", "w+");
    fprintf(fp, "#n");
    
    for(int i=1; i<M+1; i++) fprintf(fp,"\tRe(v_{%i,0})\tIm(v_{%i,0})", i, i);
    
    for(int i=-M; i<M+1; i++)
    for(int j=1; j<M+1; j++)
    {
        fprintf(fp,"\tRe(v_{%i,%i})\tIm(v_{%i,%i})", i, j, i, j);
    }
    fprintf(fp, "\n");
        
    for( int i=0; i<mcmc_trials; i++)
    {
        fprintf(fp, "%i\t", i );
        for( int j=0; j<param_size; ++j)
        {
            temp_mode = chain[ i*param_size + j ];
            fprintf(fp, "%f\t%f\t", creal(temp_mode), cimag(temp_mode) );
        }
        fprintf(fp, "\n");
    }
    fclose(fp);    
}


void output_average_file_complex( int M, int mcmc_trials, int param_size, double _Complex *chain ){
    /*
     * Moving average text file.
     *
     */
    double _Complex cumulative_sums[param_size];
    for( int i=0; i<param_size; ++i )
        cumulative_sums[i] = 0.0;
        
    FILE *fp;
    fp = fopen("Testing/MCMC/data_complex.txt", "w+");
    fprintf(fp, "#n");
    
    for(int i=1; i<M+1; i++)
    {
        fprintf(fp,"\tRe(v_{%i,0})\tIm(v_{%i,0})", i, i);
    }
    
    for(int i=-M; i<M+1; i++)
    for(int j=1; j<M+1; j++)
    {
        fprintf(fp,"\tRe(v_{%i,%i})\tIm(v_{%i,%i})", i, j, i, j);
    }
    fprintf(fp, "\n");
        
    double _Complex temp_mode;
    int trials_so_far;
    for( int i=0; i<mcmc_trials; i++)
    {
        fprintf(fp, "%i\t", i );
        for( int j=0; j<param_size; ++j)
        {
        cumulative_sums[j] += chain[ i*param_size + j ];
        fprintf(fp, "%f\t%f\t",
                    creal(cumulative_sums[j])/(i+1), 
                    cimag(cumulative_sums[j])/(i+1) );
        }
        fprintf(fp, "\n");
    }
    fclose(fp);    
}

void output_report( int M, int mcmc_trials, int param_size, double _Complex *chain, double _Complex *means, double acceptance_rate, 
                    int k, int P, double obs_sigma, double diff_sigma, double obs_delta ) // TODO: pass in settings/configuration object.
{
    FILE *fp;
    fp = fopen("Testing/MCMC/report.txt", "a");
    
    int idx;
    double sigma_squared;
    
    
    for( int j=0; j<param_size; ++j)
    {
        sigma_squared = 0;
        for( int i=0; i<mcmc_trials; ++i )
        {
            idx = i*param_size + j;
            sigma_squared += cabs( chain[idx] - means[j] )*cabs( chain[idx] - means[j] );
        }
        sigma_squared /= mcmc_trials;
        fprintf(fp, "%i\t%i\t%f\t%f\t%f\t%f\n", k, P, obs_sigma, diff_sigma, obs_delta, sigma_squared);
    }
    
    fclose(fp);
}