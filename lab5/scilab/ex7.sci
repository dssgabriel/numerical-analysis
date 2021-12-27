exec ./scilab/lib_iterative.sci

[file, mode] = mopen("data/richardson.dat", "wb");

// Richardson test
n = 64;
T0 = -5.0;
T1 = 5.0;
e = 10^(-15);
max_iter = 10000;

A = set_GB_operator_poisson1D(n);
b = set_RHS_poisson1D(n, T0, T1);
for alpha = 0.1 : 0.01 : 0.52
    [iter_r, res_r, x_r] = richardson(A, b, alpha, e, max_iter);
    mfprintf(file, "%e\t%d\t%e\n", alpha, iter_r($), res_r($));
end

mclose(file);
