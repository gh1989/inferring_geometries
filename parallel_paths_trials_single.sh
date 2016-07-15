#!/bin/bash   
# Different values of parallel_paths  
mkdir Testing/MCMC/saved_data/parallel_paths/     
for i in `seq 18 2 48`; do
  bin/mcmc_langevin_single_demo.exe -k $i -N 10000 &
  wait
  mkdir Testing/MCMC/saved_data/parallel_paths/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/parallel_paths/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/parallel_paths/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/parallel_paths/$i/data_complex_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/parallel_paths/report.txt
