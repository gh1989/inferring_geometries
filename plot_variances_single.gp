set term png
set ylabel "variance"

set xlabel "ind. paths"
set output 'Testing/MCMC/saved_data/parallel_paths/variance_plot.png'
plot 'Testing/MCMC/saved_data/parallel_paths/report.txt' using 1:6 title "Variance of Posterior Distribution"

set xlabel "path length"
set output 'Testing/MCMC/saved_data/path_length/variance_plot.png'
plot 'Testing/MCMC/saved_data/path_length/report.txt' using 2:6 title "Variance of Posterior Distribution"

set xlabel "obs. sigma"
set output 'Testing/MCMC/saved_data/obs_sigma/variance_plot.png'
plot 'Testing/MCMC/saved_data/obs_sigma/report.txt' using 3:6 title "Variance of Posterior Distribution"

set xlabel "diff. sigma"
set output 'Testing/MCMC/saved_data/diff_sigma/variance_plot.png'
plot 'Testing/MCMC/saved_data/diff_sigma/report.txt' using 4:6 title "Variance of Posterior Distribution"

set xlabel "obs. delta"
set output 'Testing/MCMC/saved_data/obs_delta/variance_plot.png'
plot 'Testing/MCMC/saved_data/obs_delta/report.txt' using 5:6 title "Variance of Posterior Distribution"

