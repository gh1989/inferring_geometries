#pragma once

#ifndef M_PI
    #define M_PI 3.14159265359
#endif

#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

using namespace Eigen;

void print_tensor( Tensor<double, 4>, int, int, int );
void print_tensor( Tensor<double, 2>, int, int );
void print_tensor( Tensor<double, 3>, int, int );