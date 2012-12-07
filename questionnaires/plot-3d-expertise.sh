#!/bin/sh
data=$1
pdf=$2
gnuplot <<EOF 
reset
set key off
set border 1
set xrange [-0.5:]
set yrange [-0.5:]
set ytics ("never heard of it" 0, "know what it is / might have used it occasionally" 1, "use it but don't really understand it" 2, "use it regularly and feel I understand it well" 3, "expert" 4)
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
