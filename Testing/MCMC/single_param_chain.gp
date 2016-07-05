set term png  
set xlabel "Timestep n"
set ylabel ""
set grid
set key outside

set output "single_param_chain_average.png"
set title "Moving Average for C"
plot "single_param_data.txt" u 1:2 with lines lc rgb "red" title "<C>", \
	"single_param_data_ts.txt" u 1:2 pt 2 lc rgb "red" title "C"