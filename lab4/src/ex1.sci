exec ./src/myldlt.sci;
exec ./src/mylu1b.sci;

step_size = 8;
max_size = 512;
nb_runs = 8;
errs_ldlt = zeros(max_size / step_size, nb_runs);
errs_lu1b = zeros(max_size / step_size, nb_runs);
errs_lusl = zeros(max_size / step_size, nb_runs);
elapsed_ldlt = zeros(max_size / step_size, nb_runs);
elapsed_lu1b = zeros(max_size / step_size, nb_runs);
elapsed_lusl = zeros(max_size / step_size, nb_runs);

[file1, mode1] = mopen("data/ldlt.dat", "wb");
[file2, mode2] = mopen("data/lu1b.dat", "wb");
[file3, mode3] = mopen("data/lusl.dat", "wb");

for i = 8 : step_size : max_size
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
        [L1, U] = mylu1b(A);
        elapsed_lu1b(i / step_size, j) = toc();
        
        // Scilab LU
        tic();
        [L3, U] = lu(A);
        elapsed_lusl(i / step_size, j) = toc();

        // Extract A from LDLT
        L2 = tril(LDLT);
        D = zeros(i, i);
        d = diag(LDLT);
        for k = 1 : i
            L2(k, k) = 1;
            D(k, k) = d(k);
        end
        LDLT = L2 * D * L2';
        
        // Extract A from LU
        LU = L1 * U;

        errs_ldlt(i / step_size, j) = norm(A - LDLT) / norm(A);
        errs_lu1b(i / step_size, j) = norm(A - LU) / norm(A);
        errs_lusl(i / step_size, j) = norm(A - L3 * U) / norm(A);
    end
    mfprintf(file1, "%d\t%e\t%e\n", i, mean(errs_ldlt, 'c')(i / step_size), mean(elapsed_ldlt, 'c')(i / step_size));
    mfprintf(file2, "%d\t%e\t%e\n", i, mean(errs_lu1b, 'c')(i / step_size), mean(elapsed_lu1b, 'c')(i / step_size));
    mfprintf(file3, "%d\t%e\t%e\n", i, mean(errs_lusl, 'c')(i / step_size), mean(elapsed_lusl, 'c')(i / step_size));
end

mclose(file1);
mclose(file2);
mclose(file3);
