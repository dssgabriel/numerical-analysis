exec ./usolve.sci;

function [x] = gausskij3b(A, b)
    n = size(A)(1);
    x = zeros(n, 1);
    
    for k = 1 : n - 1
        for i = k + 1 : n
            m_ik = A(i, k) / A(k, k);
            b(i) = b(i) - m_ik * b(k);
            for j = k + 1 : n
                A(i, j) = A(i, j) - m_ik * A(k, j);
            end
        end
    end
    
    x = usolve(A, b);
endfunction

sizes = [2; 3; 4; 8; 12; 16; 24; 32];
n = size(sizes)(1);
forward_errs = zeros(n, 1);
backward_errs = zeros(n, 1);

for i = 1 : n
    A = rand(sizes(i), sizes(i)) + ones(sizes(i), sizes(i));
    b = rand(sizes(i), 1) + ones(sizes(i), 1);

    x1 = gausskij3b(A, b);
    x2 = A \ b;
    
    forward_errs(i) = norm(x2 - x1) / norm(x2);
    backward_errs(i) = norm(b - A * x1) / (norm(A) * norm(x1));
end
