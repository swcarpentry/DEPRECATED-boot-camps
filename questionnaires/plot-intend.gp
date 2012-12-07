reset
set key off
set border 1
set xrange [-0.5:]
set yrange [-0.5:]
set ytics ("strongly agree", "agree", "neither agree or disagree", "disagree", "strongly disagree")
# "plot" data from file "using" data
#   1 means Column 1 provides X value
#   (N): means use N as Y value
#   (0.15*$N) means use 0.15 * value of column N as magnitude
datafile="data/intend.dat"
blob=0.3
plot datafile u 1:(0):(blob*$3):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(1):(blob*$4):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(2):(blob*$5):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(3):(blob*$6):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(4):(blob*$7):xticlabels(2) lc rgb "blue" pt 7 ps variable
set terminal pdf colour enhanced font "Helvetica,8" size 12cm,10cm
set output "pdfs/intend.pdf"
replot
