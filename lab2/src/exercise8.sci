function [C] = matmat3b(A, B)
    C = zeros(size(A)(2), size(B)(1));
    for i = 1 : size(A)(1)
        for j = 1 : size(B)(1)
            for k = 1 : size(A)(2)
                C(i, j) = A(i, k) * B(k, j) + C(i, j);
            end
        end
    end
endfunction

function [C] = matmat2b(A, B)
    C = zeros(size(A)(2), size(B)(1));
    for i = 1 : size(A)(1)
        for j = 1 : size(B)(2)
            C(i, j) = A(i, :) * B(:, j) + C(i, j);
        end
    end
endfunction

function [C] = matmat1b(A, B)
    C = zeros(size(A)(2), size(B)(1));
    for i = 1 : size(A)(1)
        C(i, :) = A(i, :) * B(:, :) + C(i, :);
    end
endfunction

sizes = [10; 20; 40; 80; 160; 320];
n = size(sizes)(1);
elapsed_3b = zeros(n, 1);
elapsed_2b = zeros(n, 1);
elapsed_1b = zeros(n, 1);

for i = 1 : n
    A = rand(sizes(i), sizes(i));
    B = rand(sizes(i), sizes(i));
    
    tic();
    C = matmat3b(A, B);
    elapsed_3b(i) = toc();

    tic();
    D = matmat2b(A, B);
    elapsed_2b(i) = toc();

    tic();
    E = matmat1b(A, B);
    elapsed_1b(i) = toc(); 
end
