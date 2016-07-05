#include "fourier_transform.h"

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
		{
			//printf("%f+%fi", creal(out[i*M+j]), cimag(out[i*M+j]));
			printf("%f + %fi  \t", out[i*M+j][0], out[i*M+j][1]);
		}
		printf("\n");
	}
	
	fftw_free(in);
	fftw_free(out);
	
	return 0;
}