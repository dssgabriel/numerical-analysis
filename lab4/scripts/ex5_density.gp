set term pngcairo size 1080, 720
set grid
set xtics nomirror
set ytics nomirror
set x2tics
set xlabel "Matrix density"
set x2label "Number of non-zero elements"
set ylabel "Time (s)"
set logscale y
set output "plots/ex5_density.png"
plot "data/dspmv_density.dat" u 1:3 w lp t "mydspmv time", "data/dspmv_density.dat" u 2:4 w lp t "Scilab dspmv time"
