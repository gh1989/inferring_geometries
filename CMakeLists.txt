cmake_minimum_required (VERSION 3.0)

project (InferringGeometries)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(fftw3 STATIC IMPORTED)
set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION /usr/local/lib/libfftw3.a)
add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION /usr/local/lib/libgsl.a)
add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION /usr/local/lib/libgslcblas.a)

add_library(euler_maruyama src/euler_maruyama.cxx)


add_library(fourier_series src/fourier_series.cxx)

add_library(mcmc src/mcmc.cxx)
add_library(grid src/grid.cxx)
target_link_libraries(grid fftw3)
target_link_libraries( euler_maruyama fftw3 grid )
include_directories ( include /usr/local/include/eigen3 )

add_executable( grid_demo src/grid_demo.cxx )
target_link_libraries( grid_demo fftw3 grid fourier_series )

add_executable( mcmc_ou_demo src/mcmc_ou_demo.cxx )
target_link_libraries( mcmc_ou_demo mcmc euler_maruyama fourier_series gsl gslcblas)

add_executable(langevin_demo src/langevin_demo.cxx)
target_link_libraries(langevin_demo euler_maruyama fourier_series gsl gslcblas)
target_compile_features(langevin_demo PRIVATE cxx_range_for)

add_executable(fourier_series_demo src/fourier_series_demo.cxx)
target_link_libraries(fourier_series_demo euler_maruyama fourier_series gsl gslcblas)

add_executable(mcmc_langevin_demo src/mcmc_langevin_demo.cxx)
target_link_libraries(mcmc_langevin_demo mcmc euler_maruyama fourier_series gsl gslcblas)

add_executable(mcmc_langevin_single_demo src/mcmc_langevin_single_demo.cxx)
target_link_libraries(mcmc_langevin_single_demo mcmc euler_maruyama fourier_series gsl gslcblas)

add_executable( fourier_transform_demo src/fourier_transform_demo.cxx )
target_link_libraries( fourier_transform_demo fftw3 fourier_series gsl gslcblas )
