function [L, U] = mylu1b(A)
    n = size(A)(1);
        
    for k = 1 : n - 1
        A((k + 1) : n, k) = A((k + 1) : n, k) / A(k, k);
        A((k + 1) : n, (k + 1) : n) = A((k + 1) : n, (k + 1) : n) - A((k + 1) : n, k) * A(k, (k + 1) : n);
    end
    
    U = triu(A);
    L = tril(A);
    for i = 1 : n
        L(i, i) = 1;
    end
endfunction

sizes = [2; 3; 4; 8; 12; 16; 24; 32];
len = size(sizes)(1);
nb_runs = 200;
errs = zeros(len, nb_runs);
elapsed = zeros(len, nb_runs);

[file, mode] = mopen("data/mylu1b.dat", "wb");

for i = 1 : len
    for j = 1 : nb_runs
        A = rand(sizes(i), sizes(i)) + ones(sizes(i), sizes(i));

        tic();
        [L, U] = mylu1b(A);
        elapsed(i, j) = toc();
    
        errs(i, j) = norm(A - L * U);
    end
    
    mfprintf(file, "%d\t%e\t%e\n", sizes(i), mean(errs, 'c')(i), mean(elapsed, 'c')(i));
end

mclose(file);
