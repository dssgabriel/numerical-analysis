/**
* File: tp2_poisson1D_direct.c
*
* This file contains the main function
* to solve the Poisson 1D problem.
**/

#include <stdio.h>
#include <stdlib.h>

#include "lib_poisson1D.h"

int main(int argc, char **argv)
{
    int NRHS = 1;
    int nbpoints = 102;
    int la = nbpoints - 2;
    double T0 = -5.0, T1 = 5.0;

    printf("-- Poisson 1D --\n\n");
    double *RHS = malloc(la * sizeof(double));
    double *EX_SOL = malloc(la * sizeof(double));
    double *X = malloc(la *sizeof(double));

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "RHS.dat");
    write_vec(EX_SOL, &la, "EX_SOL.dat");
    write_vec(X, &la, "X_grid.dat");

    int kv = 1, ku = 1, kl = 1;
    int lab = kv + kl + ku + 1;

    double *AB = malloc(lab * la * sizeof(double));

    // Working array for pivot used by LU Factorization
    int *ipiv = calloc(la, sizeof(int));

    int row = 0;
    int info = 0;
    if (row) { // LAPACK_ROW_MAJOR
        set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
        // write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku,
                             NRHS, AB, la, ipiv, RHS, NRHS);

    } else { // LAPACK_COL_MAJOR
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
        info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku,
                             NRHS, AB, lab, ipiv, RHS, la);
    }    
  
    printf("INFO DGBSV = %d\n", info);
    write_xy(RHS, X, &la, "SOL.dat");

    // Relative residual
    double tmp = cblas_ddot(la, RHS, 1, RHS,1);
    tmp = sqrt(tmp);
    cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
    double relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
    relres = sqrt(relres);
    relres = relres / tmp;
    printf("The relative residual error is relres = %e\n", relres);

    // Free arrays
    free(RHS);
    free(EX_SOL);
    free(X);
    free(AB);
    free(ipiv);
    printf("\n-- End --\n");

    return 0;
}
