#!/usr/bin/gnuplot

 unset xtics
set key bottom
plot 'functions.dat' u 1:2 w l lw 3 t 'Density','functions.dat' u 1:3 w l lw 3 t 'Potential','functions.dat' u 1:4 w l lw 3 t 'Result'
#    EOF
