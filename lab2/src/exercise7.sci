sizes = [10; 20; 40; 80; 160; 320];
n = size(sizes)(1);
forward_errs = zeros(n, 1);
backward_errs = zeros(n, 1);

for i = 1 : n
    A = rand(sizes(i), sizes(i));
    xex = rand(sizes(i), 1);

    b = A * xex;
    x = A \ b;

    forward_errs(i) = norm(xex - x) / norm(xex);
    backward_errs(i) = norm(b - A * x) / (norm(A) * norm(x));
end
