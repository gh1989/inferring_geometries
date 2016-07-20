#include "util.h"

/*
 * Printing helper functions.
 */
 
void print_tensor( Tensor<double,4> x, int K, int T, int N )
{

    for( int k=0; k<K; ++k )
    {
    printf("Independent path: %i... \n", k );
    for( int i=0; i<T; ++i ){
        for( int j=0; j<N; ++j )
            printf("(%.2f,%.2f)\t", x(k,i,j,0), x(k,i,j,1) );
    printf("\n");
    }
    }
}

void print_tensor( Tensor<double,2> w, int T, int N )
{
    for( int i=0; i<T; ++i ){
        for( int j=0; j<N; ++j ) printf("%.10f\t", w(i,j) );
    printf("\n");
    }
}

void print_tensor( Tensor<double,3> y, int K, int T )
{
    for( int k=0; k<K; ++k )
    {
    printf("Independent path: %i... \n ", k);
    for( int i=0; i<T; ++i )
        printf("(%.2f,%.2f)\n", y(k,i,0), y(k,i,1) );
    }
}
