#!/usr/bin/gnuplot
set datafile separator ","
set terminal png
set output "count.png"
plot "sum.csv" using 1:4 with lines title 'cellcount', "sum.csv" using 1:5 with lines title 'activated cells'
set output "eps.png"
plot "sum.csv" using 1:6 with lines title 'eps particles'
set output "cs.png"
plot "sum.csv" using 1:2 with lines title 'substrate'
set output "cq.png"
plot "sum.csv" using 1:3 with lines title 'QS molecule'
