function [A] = myldlt(A)
    [n, _] = size(A);
    piv = zeros(n, 1);

    for j = 1 : n;
        // First iteration must be done manually
        if j == 1
            A(2 : n, 1) = A(2 : n, 1) / A(1, 1);
        else
            for i = 1 : j - 1;
                piv(i) = A(j, i) * A(i, i);
            end
            A(j, j) = A(j, j) - A(j, 1 : j - 1) * piv(1 : j - 1);
            A(j + 1 : n, j) = (A(j + 1 : n, j) - A(j + 1 : n, 1 : j - 1) * piv(1 : j - 1)) / A(j, j); 
        end
    end
endfunction
