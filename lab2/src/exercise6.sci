x = [1, 2, 3, 4]
y = [1; 2; 3; 4]

// z = x + y -- Inconsistent row/column dimensions
s = x * y

x_size = size(x)
y_size = size(y)

x_norm = norm(x)
y_norm = norm(y)

A = [1, 2, 3;
     4, 5, 6;
     7, 8, 9;
     10, 11, 12]

A_transpose = A'

A = [1, 2, 3;
     4, 5, 6;
     7, 8, 9]
B = A'

C = A + B
D = A * B

A_cond = cond(A)
