set term pngcairo size 1080, 720
set grid
set xtics nomirror
set ytics nomirror
set x2tics
set y2tics
set xlabel "Matrix size"
set x2label "Number of non-zero elements"
set ylabel "Distance to solution"
set y2label "Time (s)"
set logscale x
set logscale y2
set output "plots/ex5_size.png"
plot "data/dspmv_size.dat" u 2:3 w lp t "mydspmv error" axes x1y1, "data/dspmv_size.dat" u 1:4 w lp t "mydspmv time" axes x1y2, "data/dspmv_size.dat" u 1:5 w lp t "Scilab dspmv time" axes x1y2
