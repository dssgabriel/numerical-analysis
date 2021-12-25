set term pngcairo size 1080, 720
set grid
set xlabel "Matrix density"
set ylabel "Time (s)"
# set logscale x
# set logscale y
set output "plots/ex5_density.png"
plot "data/dspmv_density.dat" u 1:2 w lp t "mydspmv time", "data/dspmv_density.dat" u 1:3 w lp t "Scilab dspmv time"
