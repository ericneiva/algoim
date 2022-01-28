set term png
set output "h_conv_surface_ellipse.png"

set title 'ellipse perimeter'
set key off
set logscale x 2
set logscale y 10
set grid

plot 'h_convergence_surface.dat' using 1:2 with points linestyle 1 title 'k = 1' ,\
     'h_convergence_surface.dat' using 1:4 with points linestyle 2 title 'k = 2' ,\
     'h_convergence_surface.dat' using 1:6 with points linestyle 3 title 'k = 4' ,\
     'h_convergence_surface.dat' using 1:8 with points linestyle 4 title 'k = 8' ,\
     'h_conv_th_surf_ellipse.dat' using 1:2 with lines linestyle 1, \
     'h_conv_th_surf_ellipse.dat' using 1:3 with lines linestyle 2, \
     'h_conv_th_surf_ellipse.dat' using 1:4 with lines linestyle 3, \
     'h_conv_th_surf_ellipse.dat' using 1:5 with lines linestyle 4