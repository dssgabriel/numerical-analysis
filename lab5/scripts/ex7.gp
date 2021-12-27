set term pngcairo size 1080, 720
set grid
set view 80, 165, 1, 1
set xlabel "alpha"
set ylabel "Number of iterations"
set zlabel "Relative residual"
set logscale z
set output "plots/ex7.png"
splot "data/richardson.dat" matrix u 2:1:3 every 1::1 w l palette t "Convergence depending on alpha"
