set term gif animate delay 8 size 600, 600
outfile = 'burn_animation.gif'
set output outfile
set view map
do for [k=0:127] {
infile = sprintf("burn_frames/burn_distribution_%04u.txt", k)
splot infile matrix with image
}
