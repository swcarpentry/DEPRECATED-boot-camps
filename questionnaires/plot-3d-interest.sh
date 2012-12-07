#!/bin/sh
data=$1
pdf=$2
gnuplot <<EOF 
reset
set key off
set border 1
set xrange [-0.5:]
set yrange [-0.5:]
set ytics ("don't know and don't care" 0, "not sure and want more" 1, "definitely" 2, "know and want more" 3, "know but don't care" 4)
file="$data"
blob=0.3
plot file u 1:(0):(blob*\$3):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     file u 1:(1):(blob*\$4):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     file u 1:(2):(blob*\$5):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     file u 1:(3):(blob*\$6):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     file u 1:(4):(blob*\$7):xticlabels(2) lc rgb "blue" pt 7 ps variable
set terminal pdf colour enhanced font "Helvetica,8" size 12cm,10cm
set output "$pdf"
replot
EOF
