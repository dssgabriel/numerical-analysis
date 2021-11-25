function [x] = usolve(U, b)
    n = size(b)(1);
    x = zeros(n, 1);
    
    x(n) = b(n) / U(n, n);
    for i = n - 1 : -1 : 1
        x(i) = (b(i) - U(i, (i + 1) : n) * x((i + 1) : n)) / U(i, i);
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
    U = triu(A);

    x1 = usolve(U, b);
    x2 = U \ b;
    
    forward_errs(i) = norm(x2 - x1) / norm(x2);
    backward_errs(i) = norm(b - A * x1) / (norm(A) * norm(x1));
end
