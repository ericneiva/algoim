set term png
set output "h_conv_ellipsoid.png"

set title 'ellipsoid volume'
set key off
set logscale x 2
set logscale y 10
set grid

plot 'h_convergence.dat' using 1:3 with points linestyle 1 title 'k = 1' ,\
     'h_convergence.dat' using 1:5 with points linestyle 2 title 'k = 2' ,\
     'h_convergence.dat' using 1:7 with points linestyle 3 title 'k = 4' ,\
     'h_convergence.dat' using 1:9 with points linestyle 4 title 'k = 8' ,\
     'h_conv_th_ellipsoid.dat' using 1:2 with lines linestyle 1, \
     'h_conv_th_ellipsoid.dat' using 1:3 with lines linestyle 2, \
     'h_conv_th_ellipsoid.dat' using 1:4 with lines linestyle 3, \
     'h_conv_th_ellipsoid.dat' using 1:5 with lines linestyle 4