#!/bin/sh
data=$1
pdf=$2
gnuplot <<EOF 
reset
set key off
set border 1
set boxwidth 1 relative
set style fill solid
set style data histogram
file="$data"
plot file using 2:xticlabels(1) notitle lc rgb "blue"
set terminal pdf colour enhanced font "Helvetica,8" size 8cm,10cm
set output "$pdf"
replot
EOF
