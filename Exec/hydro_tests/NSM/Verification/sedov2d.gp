# plot the Sod results.  Here we assume that we have a file called sedov2d.out
# with columns containing x, density, velocity, and pressure

set term post eps color solid enhanced
set output 'sedov_cyl.eps'

set multiplot;

set size 0.5, 0.5;

set xrange [0:0.5]
set xlabel "r";

set origin 0.0, 0.5;
set ylabel "density";
plot 'sedov_2d_cyl.out' using 1:2 notitle,\
     'cylindrical_sedov.dat' using 2:3 notitle with lines;

set origin 0.5, 0.5;
set ylabel "velocity";
plot 'sedov_2d_cyl.out' using 1:3 notitle,\
     'cylindrical_sedov.dat' using 2:6 notitle with lines;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'sedov_2d_cyl.out' using 1:4 notitle,\
     'cylindrical_sedov.dat' using 2:5 notitle with lines;


set origin 0.5, 0.0;
set logscale y
set format y "10^{%L}"

set ylabel "internal energy";
plot 'sedov_2d_cyl.out' using 1:5 notitle,\
     'cylindrical_sedov.dat' using 2:4 notitle with lines;

unset multiplot;
set term x11;
