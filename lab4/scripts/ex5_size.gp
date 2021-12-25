set term pngcairo size 1080, 720
set grid
set xlabel "Matrix size"
set ylabel "Time (s)"
# set logscale x
# set logscale y
set output "plots/ex5_size.png"
plot "data/dspmv_size.dat" u 1:3 w lp t "mydspmv time", "data/dspmv_size.dat" u 1:4 w lp t "Scilab dspmv time"
