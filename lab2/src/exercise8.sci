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

// Sizes of the matrices
sizes = [10; 20; 40; 80; 160; 240; 320];
len = size(sizes)(1);

nb_runs = 10;

// Matrices in which we are going to store the time taken
elapsed_3b = zeros(len, nb_runs);
elapsed_2b = zeros(len, nb_runs);
elapsed_1b = zeros(len, nb_runs);

// Opening the output file
[file, mode] = mopen("data/exercise8.dat", "wb");

for i = 1 : len
    for j = 1 : nb_runs
        // Generating matrices
        A = rand(sizes(i), sizes(i));
        B = rand(sizes(i), sizes(i));
    
        tic();
        C = matmat3b(A, B);
        elapsed_3b(i, j) = toc();

        tic();
        D = matmat2b(A, B);
        elapsed_2b(i, j) = toc();

        tic();
        E = matmat1b(A, B);
        elapsed_1b(i, j) = toc();
    end
    
    // Writing the mean compute time for matrices of size `sizes(i)` to the output
    mfprintf(file, "%d\t%e\t%e\t%e\n", sizes(i), mean(elapsed_3b, 'c')(i), mean(elapsed_2b, 'c')(i), mean(elapsed_1b, 'c')(i));
end

mclose(file);
