set term pngcairo size 1080, 720
set grid
set xlabel "Number of iterations"
set ylabel "Relative residual"
set logscale y
set output "plots/ex6_conv.png"
plot "data/conv_j_20.dat" u 1:2 w l t "Jacobi (n = 20)", "data/conv_gs_20.dat" u 1:2 w l t "Gauss-Seidel (n = 20)", "data/conv_j_21.dat" u 1:2 w l t "Jacobi (n = 21)", "data/conv_gs_21.dat" u 1:2 w l t "Gauss-Seidel (n = 21)"
