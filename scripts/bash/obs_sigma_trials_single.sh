#!/bin/bash   
# different values of obs_sigma
mkdir Testing/MCMC/saved_data/obs_sigma/
for i in `seq 0.0025 0.0025 0.05`; do
  bin/mcmc_langevin_single_demo.exe -N 10000 -d $i &
  wait
  mkdir Testing/MCMC/saved_data/obs_sigma/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/obs_sigma/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/obs_sigma/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/obs_sigma/$i/data_complex_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/obs_sigma/report.txt