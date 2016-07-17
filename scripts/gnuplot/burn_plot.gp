set term png
set output "burn_plot.png"
set view map
splot "burn_distribution.txt" matrix with image
