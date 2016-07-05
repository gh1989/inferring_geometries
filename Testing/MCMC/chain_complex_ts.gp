set term png  
set output "constants/chain_ts_complex_0.png"
set title "Time series for V_{1,0}"
set xlabel "Timestep n"
set ylabel "Mean(n)"
set grid
set key outside
plot "data_complex_ts.txt" u 1:2 title "Re", \
	"data_complex_ts.txt" u 1:3 title "Im", \
	
	set term png  
set output "constants/chain_ts_complex_1.png"
set title "Time series for V_{-1,1}"
set xlabel "Timestep n"
set ylabel "Mean(n)"
set grid
set key outside
plot "data_complex_ts.txt" u 1:4 title "Re", \
	"data_complex_ts.txt" u 1:5 title "Im", \

set term png  
set output "constants/chain_ts_complex_2.png"
set title "Time series for V_{0,1}"
set xlabel "Timestep n"
set ylabel "Mean(n)"
set grid
set key outside
plot "data_complex_ts.txt" u 1:6 title "Re", \
	"data_complex_ts.txt" u 1:7 title "Im", \
	
	set term png  
set output "constants/chain_ts_complex_3.png"
set title "Time series for V_{1,1}"
set xlabel "Timestep n"
set ylabel "Mean(n)"
set grid
set key outside
plot "data_complex_ts.txt" u 1:8 title "Re", \
	"data_complex_ts.txt" u 1:9 title "Im", \
