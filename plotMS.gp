set datafile separator ","
set logscale x
set format x "10^{%L}"
set xlabel "Input Size (# elements)"
set ylabel "Time (s)"
set title "Sequential Mergesort Performance"
set grid

plot "./seqMS.csv" using 1:2 with points pointtype 7 pointsize 1 linecolor rgb "blue" notitle
