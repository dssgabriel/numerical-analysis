#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lib_lu_tridiag.h"
#include "lib_poisson1D.h"

i32 lu_tridiag(size_t ld, size_t lab, f64 *A, f64 *x, f64 *b)
{
    if (!A)
        return printf("error: in function `lu_tridiag`\n\tmatrix `A` is a null pointer\n"), -1;
    else if (!x)
        return printf("error: in function `lu_tridiag`\n\tvector `x` is a null pointer\n"), -2;
    else if (!b)
        return printf("error: in function `lu_tridiag`\n\tvector `b` is a null pointer\n"), -3;
    else if (ld < 1)
        return printf("error: in function `lu_tridiag`\n\tscalar `ld` must be greater than 0\n"), -4;
    else if (lab != 3)
        return printf("error: in function `lu_tridiag`\n\tscalar `lab` must be equal to 3 (i.e. number of columns must be 3)\n"), -5;

    // Initialization
    for (size_t i = 1; i < ld; i++) {
        A[lab * i + 1] -= A[lab * i] * A[lab * (i - 1) + 2] / A[lab * (i - 1) + 1];
        A[lab * i] /= A[lab * (i - 1) + 1];
    }

    struct timespec st, ed;
    clock_gettime(CLOCK_MONOTONIC_RAW, &st);
    // Resolve down
    x[0] = b[0];
    for (size_t i = 1; i < ld; i++)
        x[i] = b[i] - A[lab * i] * x[i - 1];

    // Resolve up
    b[ld - 1] = x[ld - 1] / A[lab * (ld - 1) + 1];
    for (ssize_t i = ld - 2; i > -1; i--)
        b[i] = (x[i] - A[lab * i + 2] * b[i + 1]) / A[lab * i + 1];
    clock_gettime(CLOCK_MONOTONIC_RAW, &ed);
    printf("LU tridiag time = %ldns\n", ed.tv_nsec - st.tv_nsec);

    return 0;
}
