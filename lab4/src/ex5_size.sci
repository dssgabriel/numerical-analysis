exec ./src/sparse.sci

start_size = 16;
step_size = start_size;
max_size = 512;
density = 0.25;
nb_runs = 10;

errs = zeros(max_size / step_size, nb_runs);
elapsed_mydspmv = zeros(max_size / step_size, nb_runs);
elapsed_sldspmv = zeros(max_size / step_size, nb_runs);

[file, mode] = mopen("data/dspmv_size.dat", "wb");

/// Measuring time complexity when varying the matrix size
for i = start_size : step_size : max_size
    for j = 1 : nb_runs
        A = sprand(i, i, density);
        A = full(A);
        [n, m] = size(A);
        [AA, JA, IA] = sp2csr(A);
        x = ones(m, 1);

        tic();
        y = mydspmv(AA, JA, IA, x, n);
        elapsed_mydspmv(i / step_size, j) = toc();
        
        tic();
        z = A * x;
        elapsed_sldspmv(i / step_size, j) = toc();

        errs(i / step_size, j) = norm(z - y) / norm(z);
    end
    mfprintf(file, "%d\t%d\t%e\t%e\t%e\n", i, size(AA, 2), mean(errs, 'c')(i / step_size), mean(elapsed_mydspmv, 'c')(i / step_size), mean(elapsed_sldspmv, 'c')(i / step_size));
end

mclose(file);
