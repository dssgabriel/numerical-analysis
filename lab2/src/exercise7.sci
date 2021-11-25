A = rand(10000, 10000)
xex = rand(10000, 1)

b = A * xex

x = A \ b

forward_err = (norm(x) - norm(xex)) / norm(xex)
backward_err = norm(b - A * x) / (norm(A) * norm(x))
