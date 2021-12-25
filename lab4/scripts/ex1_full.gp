set term pngcairo size 1080, 720
set grid
set ytics nomirror
set y2tics
set xlabel "Matrix size"
set ylabel "Relative error"
set y2label "Time (s)"
set output "plots/ex1_full.png"
plot "data/ldlt.dat" u 1:2 w lp t "myldlt error" axes x1y1, "data/ldlt.dat" u 1:3 w lp t "myldlt time" axes x1y2, "data/lu1b.dat" u 1:2 w lp t "mylu1b error" axes x1y1, "data/lu1b.dat" u 1:3 w lp t "mylu1b time" axes x1y2, "data/lusl.dat" u 1:2 w lp t "scilab lu error" axes x1y1, "data/lusl.dat" u 1:3 w lp t "scilab lu time" axes x1y2
