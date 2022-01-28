set term png
set output "q_conv_ellipse.png"

set title 'ellipse area'
set key b l
set logscale x 2
set logscale y 10
set grid

plot 'q_convergence.dat' using 1:2 with linespoints linestyle 1 title 'n = 4' ,\
     'q_convergence.dat' using 1:4 with linespoints linestyle 2 title 'n = 32' ,\
     'q_convergence.dat' using 1:6 with linespoints linestyle 3 title 'n = 128'