set term png
set output "h_conv_surface_ellipsoid_report.png"

set title 'Ellipsoid surface area error'
set key off
set xlabel 'n'
set ylabel 'error'
set logscale x 2
set logscale y 10
set grid

plot 'h_conv_surface_ellipsoid.dat' using 1:3 with points linestyle 1 title 'k = 1' ,\
     'h_conv_surface_ellipsoid.dat' using 1:5 with points linestyle 2 title 'k = 2' ,\
     'h_conv_surface_ellipsoid.dat' using 1:7 with points linestyle 3 title 'k = 4' ,\
     'h_conv_surface_ellipsoid.dat' using 1:9 with points linestyle 4 title 'k = 8' ,\
     'h_conv_th_surf_ellipsoid.dat' using 1:2 with lines linestyle 1, \
     'h_conv_th_surf_ellipsoid.dat' using 1:3 with lines linestyle 2, \
     'h_conv_th_surf_ellipsoid.dat' using 1:4 with lines linestyle 3, \
     'h_conv_th_surf_ellipsoid.dat' using 1:5 with lines linestyle 4