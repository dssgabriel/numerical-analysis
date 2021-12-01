set term pngcairo size 1080, 720
set grid
set ytics nomirror
set y2tics
set xlabel "Matrix size"
set ylabel "Forward error"
set y2label "Backward error"
set output "plots/exercise7.png"
plot "data/exercise7.dat" u 1:2 w lp t "Forward error" axes x1y1, "data/exercise7.dat" u 1:3 w lp t "Backward error" axes x1y2
