set term pngcairo size 1080, 720
set grid
set xlabel "Matrix size"
set ylabel "Time in seconds"
set output "plots/exercise8.png"
plot "data/exercise8.dat" u 1:2 w lp t "matmat3b", "data/exercise8.dat" u 1:3 w lp t "matmat2b", "data/exercise8.dat" u 1:4 w lp t "matmat1b"
