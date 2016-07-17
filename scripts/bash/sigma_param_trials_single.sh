#!/bin/bash   
#different values of sigma_param
mkdir Testing/MCMC/saved_data/sigma_param/
for i in  `seq 0.001 0.001 0.01`; do
  bin/mcmc_langevin_single_demo.exe -x $i -N 1000 -X 1 & 
  wait
  mkdir Testing/MCMC/saved_data/sigma_param/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/sigma_param/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/sigma_param/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/sigma_param/$i/data_complex_ts.txt
  mv Testing/MCMC/data_sigma_average_ts.txt Testing/MCMC/saved_data/sigma_param/$i/data_sigma_average_ts.txt
  mv Testing/MCMC/data_sigma_ts.txt Testing/MCMC/saved_data/sigma_param/$i/data_sigma_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/sigma_param/report.txt
