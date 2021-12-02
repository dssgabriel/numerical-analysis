set term pngcairo size 1080, 720
set grid
set ytics nomirror
set y2tics
set xlabel "Matrix size"
set ylabel "Forward error"
set y2label "Backward error"
set output "plots/exercise2.png"
plot "data/usolve.dat" u 1:2 w lp t "usolve forward error" axes x1y1, "data/lsolve.dat" u 1:2 w lp t "lsolve forward error" axes x1y1, "data/usolve.dat" u 1:3 w lp t "usolve backward error" axes x1y2, "data/lsolve.dat" u 1:3 w lp t "lsolve backward error" axes x1y2
