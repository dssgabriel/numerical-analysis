function [x] = usolve(U, b)
    n = size(b)(1);
    x = zeros(n, 1);
    
    x(n) = b(n) / U(n, n);
    for i = n - 1 : -1 : 1
        x(i) = (b(i) - U(i, (i + 1) : n) * x((i + 1) : n)) / U(i, i);
    end
endfunction

sizes = [10; 20; 40; 80; 160; 240; 320; 640];
len = size(sizes)(1);
nb_runs = 3;
forward_errs = zeros(len, nb_runs);
backward_errs = zeros(len, nb_runs);

[file, mode] = mopen("data/usolve.dat", "wb");

for i = 1 : len
    for j = 1 : nb_runs
        // The bigger the matrix, the more the condition number increases.
        // To avoid that, reduce the distribution of the random values
        // in the matrix.
        // Here, the distribution has been reduced to values between
        // 1 and 2.
        A = rand(sizes(i), sizes(i)) + ones(sizes(i), sizes(i));
        b = rand(sizes(i), 1) + ones(sizes(i), 1);
        U = triu(A);

        x1 = usolve(U, b);
        x2 = U \ b;
    
        forward_errs(i, j) = norm(x2 - x1) / norm(x2);
        backward_errs(i, j) = norm(b - U * x1) / (norm(U) * norm(x1));
    end
    
    mfprintf(file, "%d\t%e\t%e\n", sizes(i), mean(forward_errs, 'c')(i), mean(backward_errs, 'c')(i));
end

mclose(file);
