exec ./src/myldlt.sci;
exec ./src/mylu1b.sci;

start_size = 8;
step_size = start_size;
max_size = 512;
nb_runs = 10;
errs_ldlt = zeros(max_size / step_size, nb_runs);
errs_lu1b = zeros(max_size / step_size, nb_runs);
errs_lusl = zeros(max_size / step_size, nb_runs);
elapsed_ldlt = zeros(max_size / step_size, nb_runs);
elapsed_lu1b = zeros(max_size / step_size, nb_runs);
elapsed_lusl = zeros(max_size / step_size, nb_runs);

[file1, mode1] = mopen("data/ldlt.dat", "wb");
[file2, mode2] = mopen("data/lu1b.dat", "wb");
[file3, mode3] = mopen("data/lusl.dat", "wb");

for i = start_size : step_size : max_size
    for j = 1 : nb_runs
        A = rand(i, i) + ones(i, i);
        // Make the matrix symmetric
        A = A + A';

        // LDLT
        tic();
        LDLT = myldlt(A);
        elapsed_ldlt(i / step_size, j) = toc();

        // My LU
        tic();
        [L1, U1] = mylu1b(A);
        elapsed_lu1b(i / step_size, j) = toc();

        // Scilab LU
        tic();
        [L2, U2] = lu(A);
        elapsed_lusl(i / step_size, j) = toc();

        // Extract A from LDLT
        L = tril(LDLT);
        D = zeros(i, i);
        d = diag(LDLT);
        for k = 1 : i
            L(k, k) = 1;
            D(k, k) = d(k);
        end
        LDLT = L * D * L';        

        errs_ldlt(i / step_size, j) = norm(A - LDLT) / norm(A);
        errs_lu1b(i / step_size, j) = norm(A - L1 * U1) / norm(A);
        errs_lusl(i / step_size, j) = norm(A - L2 * U2) / norm(A);
    end
    mfprintf(file1, "%d\t%e\t%e\n", i, mean(errs_ldlt, 'c')(i / step_size), mean(elapsed_ldlt, 'c')(i / step_size));
    mfprintf(file2, "%d\t%e\t%e\n", i, mean(errs_lu1b, 'c')(i / step_size), mean(elapsed_lu1b, 'c')(i / step_size));
    mfprintf(file3, "%d\t%e\t%e\n", i, mean(errs_lusl, 'c')(i / step_size), mean(elapsed_lusl, 'c')(i / step_size));
end

mclose(file1);
mclose(file2);
mclose(file3);
