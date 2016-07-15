#!/bin/bash   
#different values of diff_sigma
mkdir Testing/MCMC/saved_data/diff_sigma/
for i in  `seq 0.00025 0.00025 0.005`; do
  bin/mcmc_langevin_single_demo.exe -o $i -N 10000 &
  wait
  mkdir Testing/MCMC/saved_data/diff_sigma/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/diff_sigma/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/diff_sigma/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/diff_sigma/$i/data_complex_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/diff_sigma/report.txt
