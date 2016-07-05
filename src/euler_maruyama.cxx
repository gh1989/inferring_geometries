#include "euler_maruyama.h"

const int M = 1;

int main()
{
	int n = 100;

	double dt = 0.1;

	Matrix2d C;
	C << 1, 0, 0, 1;
	
	FourierSeries V(M);
	double _Complex only_mode = 1.0/(2.0*_Complex_I);
	V.set_mode( 1,1, only_mode );
	
	Vector2d *data;
	data = (Vector2d*) malloc(n*sizeof(Vector2d));

	int seed = 1234;
	gsl_rng *r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);
	
	em_overdamped_langevin( n, dt, C, &V, data, r, Vector2d(0,0) );
	free(data);
}