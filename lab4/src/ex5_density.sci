exec ./src/sparse.sci

matrix_size = 256;
start_density = 0.01;
step_density = start_density;
max_density = 0.5;
nb_runs = 10;

elapsed_mydspmv = zeros(max_density / step_density, nb_runs);
elapsed_sldspmv = zeros(max_density / step_density, nb_runs);

[file, mode] = mopen("data/dspmv_density.dat", "wb");

/// Measuring time complexity when varying density
for i = start_density : step_density : max_density
    for j = 1 : nb_runs
        A = sprand(matrix_size, matrix_size, i);
        A = full(A);
        [AA, JA, IA] = sp2csr(A);
        x = ones(matrix_size, 1);

        tic();
        y = mydspmv(AA, JA, IA, x, matrix_size);
        elapsed_mydspmv(i / step_density, j) = toc();
        
        tic();
        z = A * x;
        elapsed_sldspmv(i / step_density, j) = toc();
    end
    mfprintf(file, "%f\t%e\t%e\n", i, mean(elapsed_mydspmv, 'c')(i / step_density), mean(elapsed_sldspmv, 'c')(i / step_density));
end

mclose(file);
