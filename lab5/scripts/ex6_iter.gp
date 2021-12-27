set term pngcairo size 1080, 720
set grid
set xlabel "Poisson matrix size"
set ylabel "Number of iterations"
set output "plots/ex6_iter.png"
plot "data/jacobi.dat" u 1:2 w l t "Jacobi", "data/gauss_seidel.dat" u 1:2 w l t "Gauss-Seidel"
