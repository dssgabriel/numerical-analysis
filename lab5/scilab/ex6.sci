exec ./scilab/lib_iterative.sci

T0 = -5.0;
T1 = 5.0;
e = 10^(-12);
max_iter = 100000;

start_size = 4;
max_size = 96;

[file1, mode1] = mopen("data/jacobi.dat", "wb");
[file2, mode2] = mopen("data/gauss_seidel.dat", "wb");

[file3, mode3] = mopen("data/conv_j_20.dat", "wb");
[file4, mode4] = mopen("data/conv_gs_20.dat", "wb");

[file5, mode5] = mopen("data/conv_j_21.dat", "wb");
[file6, mode6] = mopen("data/conv_gs_21.dat", "wb");

// Iteration test
for n = start_size : max_size;
    A = set_GB_operator_poisson1D(n);
    b = set_RHS_poisson1D(n, T0, T1);

    [iter_jac, res_jac, x_jac] = jacobi(A, b, e, max_iter);
    [iter_gs, res_gs, x_gs] = gauss_seidel(A, b, e, max_iter);

    mfprintf(file1, "%d\t%e\t%e\n", n, iter_jac($), res_jac($));
    mfprintf(file2, "%d\t%e\t%e\n", n, iter_gs($), res_gs($));
end

// Convergence tests
n = 20;
A = set_GB_operator_poisson1D(n);
b = set_RHS_poisson1D(n, T0, T1);
[iter_jac, res_jac, x_jac] = jacobi(A, b, e, max_iter);
[iter_gs, res_gs, x_gs] = gauss_seidel(A, b, e, max_iter);
for i = 1 : size(iter_jac, 1)
    mfprintf(file3, "%d\t%e\n", iter_jac(i), res_jac(i));
end
for i = 1 : size(iter_gs, 1)
    mfprintf(file4, "%d\t%e\n", iter_gs(i), res_gs(i));
end

n = 21;
A = set_GB_operator_poisson1D(n);
b = set_RHS_poisson1D(n, T0, T1);
[iter_jac, res_jac, x_jac] = jacobi(A, b, e, max_iter);
[iter_gs, res_gs, x_gs] = gauss_seidel(A, b, e, max_iter);
for i = 1 : size(iter_jac, 1)
    mfprintf(file5, "%d\t%e\n", iter_jac(i), res_jac(i));
end
for i = 1 : size(iter_gs, 1)
    mfprintf(file6, "%d\t%e\n", iter_gs(i), res_gs(i));
end

mclose(file1);
mclose(file2);
mclose(file3);
mclose(file4);
mclose(file5);
mclose(file6);
