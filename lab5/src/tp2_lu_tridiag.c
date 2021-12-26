#include <stdio.h>
#include <stdlib.h>

#include "lib_poisson1D.h"
#include "lib_lu_tridiag.h"

int main()
{
    // Constants declaration
    int la = 10;
    int lab = 3;
    int kv = 0;
    f64 T0 = -5.0f, T1 = 5.0f;

    // Matrix/vector allocation
    f64 *A = malloc(la * lab * sizeof(f64));
    f64 *x = malloc(la * sizeof(f64));
    f64 *b = malloc(la * sizeof(f64));
    f64 *EX_SOL = malloc(la * sizeof(f64));

    // Setting up
    set_grid_points_1D(x, &la);
    set_dense_RHS_DBC_1D(b, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, x, &la, &T0, &T1);
    set_GB_operator_colMajor_poisson1D(A, &lab, &la, &kv);
    write_vec(x, &la, "GRID.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(b, &la, "B.dat");

    // LU tridiag
    printf("\n--- LU Tridiagonal Factorization ---\n");
    i32 info = lu_tridiag((size_t)la, (size_t)lab, A, x, b);
    write_xy(b, x, &la, "SOL.dat"); 
    printf("Info dtridiaglu = %d\n", info); 

    // Compute relative error
    // Compute norm 2 of b
    f64 tmp = cblas_ddot(la, b, 1, b, 1);
    tmp = sqrtl(tmp);
    // Compute EX_SOL - b
    cblas_daxpy(la, -1.0f, b, 1, EX_SOL, 1);
    // Compute norm 2 of EX_SOL
    f64 rel_err = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
    rel_err = sqrtl(rel_err);
    // Compute relative error
    rel_err /= tmp;
    printf("Relative error for `lu_tridiag` = %e\n", rel_err);

    printf("--- End ---\n");

    free(EX_SOL);
    free(b);
    free(x);
    free(A);

    return 0;
}
