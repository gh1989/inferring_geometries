#pragma once
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <Eigen/Dense>
#include <math.h>

#include "fourier_series.h"
#include "grid.h"


using namespace Eigen;
using Eigen::Matrix2d;

void em_ornstein_uhlenbeck( int n, double dt, double c, Vector2d *path_generated, gsl_rng *r );
int em_overdamped_langevin(int n, double dt, Matrix2d C, FourierSeries *V, Vector2d *data, gsl_rng *r, Vector2d start_from);
int em_fast_overdamped_langevin( int n, double dt, Matrix2d C, Grid *V, Vector2d *data, gsl_rng *r, Vector2d start_from);