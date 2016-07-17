#pragma once

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

#include "smc.h"
#include "fourier_series.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void particle_mcmc();