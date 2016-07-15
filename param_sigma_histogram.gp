set term png

min = 0.001
max = 0.010
step = 0.001

do for [k=0:10] {
outfile = sprintf('Testing/MCMC/saved_data/sigma_param/%.3f/sigma_hist_plot.png',min+k*step )
set output outfile
infile = sprintf('Testing/MCMC/saved_data/sigma_param/%.3f/data_sigma_ts.txt', min+k*step )
reset
n=32 #number of intervals
max=0.1 #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"
#count and plot
plot infile u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green" notitle
}

set term gif animate delay 8 size 600, 600

min_k = 0.001
max_k = 0.010
step = 0.001
outfile = 'Testing/MCMC/saved_data/sigma_param/sigma_hist_plot.gif'


reset
n=32 #number of intervals
max=0.1 #max value
min=0. #min value
width=(max-min)/n #interval width
#function used to map a value to the intervals
hist(x,width)=width*floor(x/width)+width/2.0
set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
#set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
set boxwidth width*0.9
set style fill solid 0.5 #fillstyle
set tics out nomirror
set xlabel "x"
set ylabel "Frequency"
set output outfile

do for [k=0:9] {
infile = sprintf('Testing/MCMC/saved_data/sigma_param/%.3f/data_sigma_ts.txt', min_k+k*step )
#count and plot
set title sprintf("param sigma: %f", min_k+k*step )
plot infile u (hist($1,width)):(1.0) smooth freq w boxes lc rgb"green"
}