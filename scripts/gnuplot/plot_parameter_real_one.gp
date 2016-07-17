set term png  
set xlabel "Timestep n"
set ylabel ""
set grid
set key outside

set output "Testing/MCMC/chain_average.png"
set title "Moving Average for C"
plot "Testing/MCMC/data.txt" u 1:2 with lines lc rgb "red" title "<C>", \
	"Testing/MCMC/data_ts.txt" u 1:2 pt 2 lc rgb "red" title "C"