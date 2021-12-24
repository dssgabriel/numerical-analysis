exec ./src/sparse.sci

start_size = 16;
step_size = start_size;
max_size = 512;
density = 0.25;
nb_runs = 10;

errs = zeros(max_size / step_size, nb_runs);
elapsed_mydspmv = zeros(max_size / step_size, nb_runs);
elapsed_sldspmv = zeros(max_size / step_size, nb_runs);

[file1, mode1] = mopen("data/mydspmv_size.dat", "wb");
[file2, mode2] = mopen("data/sldspmv_size.dat", "wb");

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
    mfprintf(file1, "%d\t%e\t%e\n", i, mean(errs, 'c')(i / step_size), mean(elapsed_mydspmv, 'c')(i / step_size));
    mfprintf(file2, "%d\t%e\n", i, mean(elapsed_sldspmv, 'c')(i / step_size));
end

step_size = 0.02;
max_density = 0.48;
elapsed_mydspmv = zeros(max_density / step_size, nb_runs);
elapsed_sldspmv = zeros(max_density / step_size, nb_runs);
[file3, mode3] = mopen("data/mydspmv_density.dat", "wb");
[file4, mode4] = mopen("data/sldspmv_density.dat", "wb");

/// Testing time complexity when varying density
for i = 0.16 : step_size : max_density
    for j = 1 : nb_runs
        A = sprand(384, 384, i);
        A = full(A);
        [AA, JA, IA] = sp2csr(A);
        x = ones(384, 1);

        tic();
        y = mydspmv(AA, JA, IA, x, 384);
        elapsed_mydspmv(i / step_size, j) = toc();
        
        tic();
        z = A * x;
        elapsed_sldspmv(i / step_size, j) = toc();
    end
    mfprintf(file3, "%f\t%e\n", size(AA, 2), mean(elapsed_mydspmv, 'c')(i / step_size));
    mfprintf(file4, "%f\t%e\n", size(AA, 2), mean(elapsed_sldspmv, 'c')(i / step_size));
end

mclose(file1);
mclose(file2);
mclose(file3);
mclose(file4);
