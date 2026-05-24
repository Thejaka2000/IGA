
set hidden3d
set xlabel "x"
set ylabel "y"
set zlabel "z"

splot "surface.dat" using 1:2:3 with lines lw 1
