/**
* File: tp2_poisson1D_direct.c
*
* This file contains the main function
* to solve the Poisson 1D problem.
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "lib_poisson1D.h"

int main(int argc, char **argv)
{
    int NRHS = 1;
    int nbpoints = 12;
    int la = nbpoints - 2;
    double T0 = -5.0, T1 = 5.0;

    printf("-- Poisson 1D --\n\n");
    double *RHS = malloc(la * sizeof(double));
    double *EX_SOL = malloc(la * sizeof(double));
    double *X = malloc(la * sizeof(double));
    // Vectors for cblas_dgbmv
    double *Y_dgbmv = malloc(la * sizeof(double));
    double *RHS_dgbmv = malloc(la * sizeof(double));

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(RHS, &la, &T0, &T1);
    // Exact solution used for cblas_dgbmv
    set_dense_RHS_DBC_1D(RHS_dgbmv, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);

    write_vec(RHS, &la, "data/RHS.dat");
    // Output the exact solution used for cblas_dgbmv to a file
    write_vec(RHS_dgbmv, &la, "data/RHS_dgbmv.dat");
    write_vec(EX_SOL, &la, "data/EX_SOL.dat");
    write_vec(X, &la, "data/X_grid.dat");

    int kv = 1, ku = 1, kl = 1;
    int lab = kv + kl + ku + 1;

    double *AB = malloc(lab * la * sizeof(double));

    // Working array for pivot used by LU Factorization
    int *ipiv = calloc(la, sizeof(int));

    int row = 0;
    int info = 0;
    struct timespec before, after;

    if (row) { // LAPACK_ROW_MAJOR
        /// Exercise 3: dgbsv
        set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);
        // write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "data/AB_row.dat");
        clock_gettime(CLOCK_MONOTONIC_RAW, &before);
        info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku,
                             NRHS, AB, la, ipiv, RHS, NRHS);
        clock_gettime(CLOCK_MONOTONIC_RAW, &after);
    } else { // LAPACK_COL_MAJOR
        /// Exercise 3: dgbsv
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "data/AB_col.dat");
        clock_gettime(CLOCK_MONOTONIC_RAW, &before);
        info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku,
                             NRHS, AB, lab, ipiv, RHS, la);
        clock_gettime(CLOCK_MONOTONIC_RAW, &after);

        /// Exercise 4: cblas_dgbmv
        // Resetting matrix AB before calling CBLAS
        free(AB);
        // Removing value of kv from lab
        lab -= kv;
        kv = 0;
        // Reallocating with correct size
        AB = malloc(lab * la * sizeof(double));
        set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
        cblas_dgbmv(CblasColMajor,
                    CblasNoTrans, la, la,
                    kl, ku,
                    1.0f, AB, lab,   // alpha = 1.0f
                    EX_SOL, 1, 0.0f, // beta = 0.0f
                    Y_dgbmv, 1);
        write_vec(Y_dgbmv, &la, "data/Y_dgbmv_col.dat");
    }
  
    printf("INFO DGBSV = %d\n", info);
    printf("TIME DGBSV = %lfs\n", (after.tv_sec - before.tv_sec) + (after.tv_nsec - before.tv_nsec) / 1e9);
    write_xy(RHS, X, &la, "data/SOL.dat");

    // Relative residual
    double tmp = cblas_ddot(la, RHS, 1, RHS, 1);
    tmp = sqrt(tmp);
    cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
    double relres = cblas_ddot(la, EX_SOL, 1, EX_SOL, 1);
    relres = sqrt(relres);
    relres = relres / tmp;
    printf("The relative residual error is relres = %e\n", relres);

    /// Exercise 4: Relative error for cblas_dgbmv
    printf("\nCBLAS_DGBMV:\n");

    // Compute norm 2 of Y_dgbmv
    tmp = cblas_ddot(la, Y_dgbmv, 1, Y_dgbmv, 1);
    tmp = sqrt(tmp);

    // Compute B - Y_dgbmv
    cblas_daxpy(la, -1.0, Y_dgbmv, 1, RHS_dgbmv, 1);

    // Compute norm 2 of B
    relres = cblas_ddot(la, RHS_dgbmv, 1, RHS_dgbmv, 1);
    relres = sqrt(relres);

    // Compute relative error
    relres /= tmp;
    printf("Relative error for `cblas_dgbmv` = %e\n", relres);

    // Free arrays
    free(RHS);
    free(RHS_dgbmv);
    free(EX_SOL);
    free(X);
    free(Y_dgbmv);
    free(AB);
    free(ipiv);
    printf("\n-- End --\n");

    return 0;
}
