set term png 
set key top
set xrange[-1:1]
set yrange[-1:1] 
set object 1 rectangle from graph 0, graph 0 to graph 1, graph 1 behind fc rgbcolor 'black' fs noborder
set pm3d at b      # draw on bottom, not as 3d surface
set view map       # don't do a 3-d looking plot
set dgrid 100,100  # grid of 100x100 pixels

outfile = 'Testing/MCMC/histos/hist_plot.png'
set output outfile
infile = 'Testing/MCMC/data_posterior_density_complex.txt'
titletext = ""
set label 11 titletext at -0.8, 0.8 front nopoint tc lt 1
splot infile u 1:2:3 w pm3d title titletext