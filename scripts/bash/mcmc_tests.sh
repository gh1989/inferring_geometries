#!/bin/bash   

# different values of parallel paths  
mkdir Testing/MCMC/saved_data/parallel_paths/     
for i in 2 4 8 16 32 64; do
  bin/mcmc_langevin_demo.exe -k $i -N 10000 &
  wait
  # gnuplot histo_plot.gp
  # wait
  mkdir Testing/MCMC/saved_data/parallel_paths/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/parallel_paths/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/parallel_paths/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/parallel_paths/$i/data_complex_ts.txt
  # mv Testing/MCMC/histos/histo_plot_1.png Testing/MCMC/saved_data/parallel_paths/$i/histo_plot_1.png
  # mv Testing/MCMC/histos/histo_plot_2.png Testing/MCMC/saved_data/parallel_paths/$i/histo_plot_2.png
  # mv Testing/MCMC/histos/histo_plot_3.png Testing/MCMC/saved_data/parallel_paths/$i/histo_plot_3.png
  # mv Testing/MCMC/histos/histo_plot_4.png Testing/MCMC/saved_data/parallel_paths/$i/histo_plot_4.png
done

# different values of diffusion coefficient
mkdir Testing/MCMC/saved_data/diffusion_coeff/
for i in 0.0001 0.001 0.01 0.1 1.0; do
  bin/mcmc_langevin_demo.exe -o $i -N 10000 &
  wait
  # gnuplot histo_plot.gp
  # wait
  mkdir Testing/MCMC/saved_data/diffusion_coeff/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/diffusion_coeff/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/diffusion_coeff/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/diffusion_coeff/$i/data_complex_ts.txt
  # mv Testing/MCMC/histos/histo_plot_1.png Testing/MCMC/saved_data/diffusion_coeff/$i/histo_plot_1.png
  # mv Testing/MCMC/histos/histo_plot_2.png Testing/MCMC/saved_data/diffusion_coeff/$i/histo_plot_2.png
  # mv Testing/MCMC/histos/histo_plot_3.png Testing/MCMC/saved_data/diffusion_coeff/$i/histo_plot_3.png
  # mv Testing/MCMC/histos/histo_plot_4.png Testing/MCMC/saved_data/diffusion_coeff/$i/histo_plot_4.png
done

# different values of path length
mkdir Testing/MCMC/saved_data/path_length/
for i in 2 4 8 16 32 64; do
  bin/mcmc_langevin_demo.exe -P $i -N 10000 &
  wait
  # gnuplot histo_plot.gp
  # wait
  mkdir Testing/MCMC/saved_data/path_length/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/path_length/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/path_length/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/path_length/$i/data_complex_ts.txt
  # mv Testing/MCMC/histos/histo_plot_1.png Testing/MCMC/saved_data/path_length/$i/histo_plot_1.png
  # mv Testing/MCMC/histos/histo_plot_2.png Testing/MCMC/saved_data/path_length/$i/histo_plot_2.png
  # mv Testing/MCMC/histos/histo_plot_3.png Testing/MCMC/saved_data/path_length/$i/histo_plot_3.png
  # mv Testing/MCMC/histos/histo_plot_4.png Testing/MCMC/saved_data/path_length/$i/histo_plot_4.png
done

# different values of obs delta
mkdir Testing/MCMC/saved_data/obs_delta/
for i in 0.0001 0.001 0.01 0.1 1.0; do
  bin/mcmc_langevin_demo.exe -N 10000 -d $i &
  wait
  # gnuplot histo_plot.gp
  # wait
  mkdir Testing/MCMC/saved_data/obs_delta/$i/
  mv Testing/MCMC/data_posterior_density_complex.txt Testing/MCMC/saved_data/obs_delta/$i/data_posterior_density_complex.txt
  mv Testing/MCMC/data_complex.txt Testing/MCMC/saved_data/obs_delta/$i/data_complex.txt
  mv Testing/MCMC/data_complex_ts.txt Testing/MCMC/saved_data/obs_delta/$i/data_complex_ts.txt
  # mv Testing/MCMC/histos/histo_plot_1.png Testing/MCMC/saved_data/obs_delta/$i/histo_plot_1.png
  # mv Testing/MCMC/histos/histo_plot_2.png Testing/MCMC/saved_data/obs_delta/$i/histo_plot_2.png
  # mv Testing/MCMC/histos/histo_plot_3.png Testing/MCMC/saved_data/obs_delta/$i/histo_plot_3.png
  # mv Testing/MCMC/histos/histo_plot_4.png Testing/MCMC/saved_data/obs_delta/$i/histo_plot_4.png
done
