set datafile separator ","
set logscale x
set format x "10^{%L}"
set xlabel "Input Size (#elements)"
set ylabel "SpeedUp"
set title "Mergesort SpeedUp Performance"
set grid
set key left

plot "./speedUp.csv" using 1:2 with linespoints pointtype 7 pointsize 1 linecolor rgb "red" notitle 
