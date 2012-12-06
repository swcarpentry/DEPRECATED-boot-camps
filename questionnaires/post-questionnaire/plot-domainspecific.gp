reset
set key off
set border 1
set boxwidth 1 relative
set style fill solid
set style data histogram
plot "data/domainspecific.dat" using 2:xticlabels(1) notitle lc rgb "blue"
set terminal pdf colour enhanced font "Helvetica,8" size 8cm,10cm
set output "pdfs/domainspecific.pdf"
replot
