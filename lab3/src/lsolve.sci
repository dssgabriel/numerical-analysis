function [x] = lsolve(L, b)
    n = size(b)(1)
    x = zeros(n, 1);
    
    x(1) = b(1) / L(1, 1);
    for i = 2 : n
        x(i) = (b(i) - L(i, 1 : (i - 1)) * x(1 : (i - 1))) / L(i, i);
    end
endfunction

sizes = [10; 20; 40; 80; 160; 320];
n = size(sizes)(1);
forward_errs = zeros(n, 1);
backward_errs = zeros(n, 1);

for i = 1 : n
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
    
    forward_errs(i) = norm(x2 - x1) / norm(x2);
    backward_errs(i) = norm(b - A * x1) / (norm(A) * norm(x1));
end
