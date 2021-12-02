function [x] = lsolve(L, b)
    n = size(b)(1)
    x = zeros(n, 1);
    
    x(1) = b(1) / L(1, 1);
    for i = 2 : n
        x(i) = (b(i) - L(i, 1 : (i - 1)) * x(1 : (i - 1))) / L(i, i);
    end
endfunction

sizes = [10; 20; 40; 80; 160; 240; 320; 640];
len = size(sizes)(1);
nb_runs = 3;
forward_errs = zeros(len, nb_runs);
backward_errs = zeros(len, nb_runs);

[file, mode] = mopen("data/lsolve.dat", "wb");

for i = 1 : len
    for j = 1 : nb_runs
        // The bigger the matrix, the more the condition number increases.
        // To avoid that, reduce the distribution of the random values
        // in the matrix.
        // Here, the distribution has been reduced to values between
        // 1 and 2.
        A = rand(sizes(i), sizes(i)) + ones(sizes(i), sizes(i));
        b = rand(sizes(i), 1) + ones(sizes(i), 1);
        L = tril(A);

        x1 = lsolve(L, b);
        x2 = L \ b;
    
        forward_errs(i, j) = norm(x2 - x1) / norm(x2);
        backward_errs(i, j) = norm(b - L * x1) / (norm(L) * norm(x1));
    end
    
    mfprintf(file, "%d\t%e\t%e\n", sizes(i), mean(forward_errs, 'c')(i), mean(backward_errs, 'c')(i));
end

mclose(file);

