#!/bin/bash   
# different values of obs_sigma
mkdir Testing/MCMC/saved_data/obs_delta/
for i in `seq 1 8`; do
  j=`python <<-END
		j = 32 * 0.001 / ( 2*$i )
		print j
	END`
  echo $j
  bin/mcmc_langevin_single_demo.exe -N 1000 -P 32 -j $j -d 0.001 &
  wait
  mkdir Testing/MCMC/saved_data/obs_delta/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/obs_delta/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/obs_delta/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/obs_delta/$i/data_complex_ts.txt
done
mv Testing/MCMC/report.txt Testing/MCMC/saved_data/obs_delta/report.txt