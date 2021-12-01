// 1. Declaring a line vector.
x = [1, 2, 3, 4];

// 2. Declaring a column vector.
y = [1; 2; 3; 4];

// 3. Basic vector operations.
//z = x + y; -- error: inconsistent row/column dimensions
s = x * y;

// 4. Computing the size of vectors.
x_size = size(x);
y_size = size(y);

// 5. Computing the norm of vectors (2 by default).
x_norm = norm(x);
y_norm = norm(y);

// 6. Declaring a matrix.
A = [1, 2, 3;
     4, 5, 6;
     7, 8, 9;
     10, 11, 12];

// 7. Computing the transpose of a matrix.
A_transpose = A';

// Declaring square matrices
A = [1, 2, 3;
     4, 5, 6;
     7, 8, 9];
B = A';

// 8. Basic matrix operations.
C = A + B;
D = A * B;

// 9. Computing condition number of a square matrix.
A_cond = cond(A);
