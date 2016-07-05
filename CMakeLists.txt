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

add_library(euler src/euler_maruyama.cxx)
add_library(fourier src/fourier_series.cxx)
add_library(mcmc src/mcmc.cxx)
add_library(grid src/grid.cxx)

include_directories ( include /usr/include/eigen3)

add_executable( grid_test src/grid.cxx )
target_link_libraries( grid_test fftw3 fourier  )

#add_executable( ou_mcmc src/ou_mcmc.cxx )
#target_link_libraries( ou_mcmc euler fourier gsl gslcblas)

#add_executable( euler_maruyama src/euler_maruyama.cxx )
#target_link_libraries( euler_maruyama euler fourier gsl gslcblas)

#add_executable(trajectory src/trajectory.cxx)
#target_link_libraries(trajectory fftw3 euler fourier gsl gslcblas)
#target_compile_features(trajectory PRIVATE cxx_range_for)

#add_executable(fourier_series src/fourier_series.cxx)
#target_link_libraries(fourier_series fftw3 euler fourier gsl gslcblas)

#add_executable(mcmcfuncs src/mcmc.cxx)
#target_link_libraries(mcmcfuncs fftw3 euler fourier gsl gslcblas)

#add_executable(one_param_mcmc src/one_param_mcmc.cxx)
#target_link_libraries(one_param_mcmc fftw3 euler fourier gsl gslcblas)

add_executable( fourier_transform src/fourier_transform.cxx )
target_link_libraries( fourier_transform fftw3 gsl gslcblas )
