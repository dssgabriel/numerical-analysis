set term pngcairo size 1080, 720
set grid
set xlabel "Matrix size"
set ylabel "Error"
set output "plots/exercise4.png"
plot "data/mylu3b.dat" w lp t "Error"
