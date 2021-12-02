function [L, U] = mylu3b(A)
    n = size(A)(1);
        
    for k = 1 : n - 1
        for i = k + 1 : n
            A(i, k) = A(i, k) / A(k, k);
        end
        for i = k + 1 : n
            for j = k + 1 : n
                A(i, j) = A(i, j) - A(i, k) * A(k, j);
            end
        end
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

[file, mode] = mopen("data/mylu3b.dat", "wb");

for i = 1 : len
    for j = 1 : nb_runs
        A = rand(sizes(i), sizes(i)) + ones(sizes(i), sizes(i));

        tic();
        [L, U] = mylu3b(A);
        elapsed(i, j) = toc();
    
        errs(i, j) = norm(A - L * U);
    end
    
    mfprintf(file, "%d\t%e\t%e\n", sizes(i), mean(errs, 'c')(i), mean(elapsed, 'c')(i));
end

mclose(file);
