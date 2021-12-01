// Sizes of the matrix
sizes = [100; 1000; 10000];
len = size(sizes)(1);

// Array in which we are going to store the errors
forward_errs = zeros(len, 1);
backward_errs = zeros(len, 1);

// Opening the file
[file, mode] = mopen("data/exercise7.dat", "wb");

for i = 1 : len
    // Generating the linear system
    A = rand(sizes(i), sizes(i));
    xex = rand(sizes(i), 1);

    b = A * xex;
    x = A \ b;

    forward_errs(i) = norm(xex - x) / norm(xex);
    backward_errs(i) = norm(b - A * x) / (norm(A) * norm(x));
    
    // Writing the computed errors for systems of size `sizes(i)` to the output
    mfprintf(file, "%d\t%e\t%e\n", sizes(i), forward_errs(i), backward_errs(i));
end

mclose(file);
