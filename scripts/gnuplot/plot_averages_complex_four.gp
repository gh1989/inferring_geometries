set term png  
set xlabel "Timestep n"
set ylabel ""
set grid
set key outside

set output "Testing/MCMC/constants/chain_average_complex_0.png"
set title "Moving Average for V_{1,0}"
plot "Testing/MCMC/data_complex.txt" u 1:2 with lines lc rgb "red" title "Re(<v>)", \
	"Testing/MCMC/data_complex.txt" u 1:3 with lines lc rgb "purple" title "Im(<v>)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:2 pt 2 lc rgb "red" title "Re(v)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:3 pt 2 lc rgb "purple" title "Im(v)", \
	

set output "Testing/MCMC/constants/chain_average_complex_1.png"
set title "Moving Average for V_{-1,1}"
plot "Testing/MCMC/data_complex.txt" u 1:4 with lines lc rgb "red" title "Re(<v>)", \
	"Testing/MCMC/data_complex.txt" u 1:5 with lines lc rgb "purple" title "Im(<v>)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:4 pt 2 lc rgb "red" title "Re(v)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:5 pt 2 lc rgb "purple" title "Im(v)", \

 
set output "Testing/MCMC/constants/chain_average_complex_2.png"
set title "Moving Average for V_{0,1}"
plot "Testing/MCMC/data_complex.txt" u 1:6 with lines lc rgb "red" title "Re(<v>)", \
	"Testing/MCMC/data_complex.txt" u 1:7 with lines lc rgb "purple" title "Im(<v>)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:6 pt 2 lc rgb "red" title "Re(v)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:7 pt 2 lc rgb "purple" title "Im(v)", \
	
 
set output "Testing/MCMC/constants/chain_average_complex_3.png"
set title "Moving Average for V_{1,1}"
plot "Testing/MCMC/data_complex.txt" u 1:8 with lines lc rgb "red" title "Re(<v>)", \
	"Testing/MCMC/data_complex.txt" u 1:9 with lines lc rgb "purple" title "Im(<v>)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:8 pt 2 lc rgb "red" title "Re(v)", \
	"Testing/MCMC/data_complex_ts.txt" u 1:9 pt 2 lc rgb "purple" title "Im(v)", \
