function [A] = set_GB_operator_poisson1D(n)
    A = 2 * eye(n, n);
    for i = 1 : n - 1
        A(i, i + 1) = -1;
        A(i + 1, i) = -1; 
    end
endfunction

function [b] = set_RHS_poisson1D(n, T0, T1)
    b = zeros(n, 1);
    b(1) = T0;
    b(n) = T1;
endfunction

function [iter, residual, x] = jacobi(A, b, e, max_iter)
    n = size(A, 1);
    x = zeros(n, 1);
    inv_D = zeros(n, n);
    for i = 1 : n
        inv_D(i, i) = 1 / A(i, i); // == 1/2 pour Poisson 1D
    end

    residual = [norm(b - A * x, 2) / norm(b)];
    counter = 1;
    iter = [counter];
    while e < residual(counter) && counter < max_iter
        x = x + inv_D * (b - A * x);
        residual = [residual; norm(b - A * x, 2) / norm(b)];
        counter = counter + 1;
        iter = [iter; counter];
    end
endfunction

function [iter, residual, x] = gauss_seidel(A, b, e, max_iter)
    n = size(A, 1);
    x = zeros(n, 1);
    D = zeros(n, n);
    for i = 1 : n
        D(i, i) = A(i, i);
    end
    E = -tril(A, -1);
    inv_DE = inv(D - E);

    residual = [norm(b - A * x, 2) / norm(b, 2)];
    counter = 1;
    iter = [counter];
    while e < residual(counter) && counter < max_iter
        x = x + inv_DE * (b - A * x);
        residual = [residual; norm(b - A * x, 2) / norm(b, 2)];
        counter = counter + 1;
        iter = [iter; counter + 1];
    end
endfunction

function [iter, residual, x] = richardson(A, b, alpha, e, max_iter)
    n = size(A, 1);
    x = zeros(n, 1);

    residual = [norm(b - A * x, 2) / norm(b)];
    counter = 1;
    iter = [counter];
    while e < residual(counter) && counter < max_iter
        x = x + alpha * (b - A * x);
        residual = [residual; norm(b - A * x, 2) / norm(b)];
        counter = counter + 1;
        iter = [iter; counter];
    end
endfunction
