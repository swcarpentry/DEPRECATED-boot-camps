reset
set key off
set border 1
set xrange [-0.5:]
set yrange [-0.5:]
set ytics ("never heard of it" 0, "know what it is / might have used it occasionally" 1, "use it but don't really understand it" 2, "use it regularly and feel I understand it well" 3, "expert" 4)

# "already know it but want to know how to use it better" 

set ytics ("don't know about it and don't care" 0, "not sure but feel I should find out more about it" 1, "definitely want to know more about it" 2, "already know it but want to ... use it better" 3, "already know about it so don't want to know more " 4)
# "plot" data from file "using" data
#   1 means Column 1 provides X value
#   (N): means use N as Y value
#   (0.15*$N) means use 0.15 * value of column N as magnitude
datafile="data/interest.dat"
blob=0.3
plot datafile u 1:(0):(blob*$3):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(1):(blob*$4):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(2):(blob*$5):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(3):(blob*$6):xticlabels(2) lc rgb "blue" pt 7 ps variable, \
     datafile u 1:(4):(blob*$7):xticlabels(2) lc rgb "blue" pt 7 ps variable
set terminal pdf colour enhanced font "Helvetica,8" size 12cm,10cm
set output "pdfs/interest.pdf"
replot
