set term pngcairo size 1080, 720
set grid
set ytics nomirror
set y2tics
set xlabel "Matrix size"
set ylabel "Error"
set y2label "Time (s)"
set output "plots/exercise4.png"
plot "data/mylu3b.dat" u 1:2 w lp t "mylu3b error" axes x1y1, "data/mylu1b.dat" u 1:2 w lp t "mylu1b error" axes x1y1, "data/mylu3b.dat" u 1:3 w lp t "mylu3b time" axes x1y2, "data/mylu1b.dat" u 1:3 w lp t "mylu1b time" axes x1y2
