function [AA, JA, IA] = sp2csr(A)
    [n, m] = size(A);
    NNZ = 1;
    AA = [];
    JA = [];
    IA = [];

    IA = [IA NNZ];
    for i = 1 : n
        for j = 1 : m
            if A(i, j) ~= 0
                AA = [AA A(i, j)];
                JA = [JA j];
                NNZ = NNZ + 1;
            end
        end
        IA = [IA NNZ];
    end
endfunction

function [y] = mydspmv(AA, JA, IA, x, nb_rows)
    y = zeros(nb_rows, 1);
    for i = 1 : nb_rows
        for j = IA(i) : IA(i + 1) - 1
            y(i) = y(i) + AA(j) * x(JA(j));
        end
    end
endfunction
