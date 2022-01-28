set term png
set output "q_conv_ellipsoid.png"

set title 'ellipsoid volume'
set key b l
set logscale x 2
set logscale y 10
set grid

plot 'q_convergence.dat' using 1:3 with linespoints linestyle 1 title 'n = 4' ,\
     'q_convergence.dat' using 1:5 with linespoints linestyle 2 title 'n = 32' ,\
     'q_convergence.dat' using 1:7 with linespoints linestyle 3 title 'n = 128'