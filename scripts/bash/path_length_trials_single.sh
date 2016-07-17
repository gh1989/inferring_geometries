#!/bin/bash   
# different values of path_length
mkdir Testing/MCMC/saved_data/path_length/
for i in `seq 2 2 32`; do
  bin/mcmc_langevin_single_demo.exe -P $i -N 10000 &
  wait
  mkdir Testing/MCMC/saved_data/path_length/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/path_length/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/path_length/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/path_length/$i/data_complex_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/path_length/report.txt
