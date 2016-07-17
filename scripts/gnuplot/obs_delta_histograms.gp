set term png 
set key top
set xrange[-1:1]
set yrange[-1:1] 
set object 1 rectangle from graph 0, graph 0 to graph 1, graph 1 behind fc rgbcolor 'black' fs noborder
set pm3d at b      # draw on bottom, not as 3d surface
set view map       # don't do a 3-d looking plot
set dgrid 100,100  # grid of 100x100 pixels

min_k = 0.001
max_k = 0.008
step = 0.001

do for [k=0:7] {
outfile = sprintf('Testing/MCMC/saved_data/obs_delta/%.3f/hist_plot.png',min_k+k*step)
set output outfile
infile = sprintf('Testing/MCMC/saved_data/obs_delta/%.3f/data_posterior_density_complex.txt',min_k+k*step)
titletext = sprintf('obs\_delta=%.3f',min_k+k*step)
set label 11 titletext at -0.8, 0.8 front nopoint tc lt 1
splot infile u 1:2:3 w pm3d title titletext
}

set term gif animate delay 8 size 600, 600
outfile = sprintf('Testing/MCMC/saved_data/obs_delta/animation_%.3f_%.3f.gif',min_k, max_k )
set output outfile

do for [k=0:7] {
infile = sprintf('Testing/MCMC/saved_data/obs_delta/%.3f/data_posterior_density_complex.txt',min_k+k*step)
titletext = sprintf('obs\_delta=%.3f',min_k+k*step)
set label 11 titletext at -0.8, 0.8 front nopoint tc lt 1
splot infile u 1:2:3 w pm3d title titletext
}
