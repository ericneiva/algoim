set term png
#set output "geo_translation_ellipsoid_volume_n=64.png"
set output "geo_translation_ellipsoid_surface_area_n=64.png"

#set title 'ellipsoid volume'
set title 'ellipsoid surface area'
set key b l
set logscale y 10
set grid

#plot 'geo_translation_n=64.out' using 1:2 with lines linestyle 1 title 'q = 4' ,\
#     'geo_translation_n=64.out' using 1:4 with lines linestyle 2 title 'q = 7' ,\
#     'geo_translation_n=64.out' using 1:6 with lines linestyle 3 title 'q = 10'
plot 'geo_translation_n=64.out' using 1:3 with lines linestyle 1 title 'q = 4' ,\
     'geo_translation_n=64.out' using 1:5 with lines linestyle 2 title 'q = 7' ,\
     'geo_translation_n=64.out' using 1:7 with lines linestyle 3 title 'q = 10'