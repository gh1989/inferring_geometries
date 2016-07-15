set term png  
set pm3d at b      # draw on bottom, not as 3d surface
set view map       # don't do a 3-d looking plot
set dgrid 100,100  # grid of 100x100 pixels

set output "Testing/MCMC/histos/histo_plot_1.png"
splot 'Testing/MCMC/data_posterior_density_complex.txt' u 1:2:3 w pm3d 
set output "Testing/MCMC/histos/histo_plot_2.png"
splot 'Testing/MCMC/data_posterior_density_complex.txt' u 4:5:6 w pm3d 
set output "Testing/MCMC/histos/histo_plot_3.png"
splot 'Testing/MCMC/data_posterior_density_complex.txt' u 7:8:9 w pm3d 
set output "Testing/MCMC/histos/histo_plot_4.png"
splot 'Testing/MCMC/data_posterior_density_complex.txt' u 10:11:12 w pm3d 