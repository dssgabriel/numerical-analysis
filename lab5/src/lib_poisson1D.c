/**
 * File: lib_poisson1D.c
 * Numerical library developed to solve 1D
 * Poisson problem (Heat equation)
 **/

#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "lib_poisson1D.h"

void set_GB_operator_rowMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
    size_t kk = (*kv) * (*la);
    if (*kv >= 0) {
        for (size_t ii = 0; ii < kk; ii++)
            AB[ii] = 0.0;
    }

    for (size_t ii = 0; ii < (*la); ii++) {
        AB[kk + ii] = -1.0;
        AB[kk + ii + (*la)] = 2.0;
        AB[kk + ii + 2 * (*la)] = -1.0;
    }

    AB[kk] = 0.0;
    AB[kk + 3 * (*la) - 1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
    for (size_t jj = 0; jj < (*la); jj++) {
        size_t kk = jj * (*lab);
        if (*kv >= 0) {
            for (size_t ii = 0; ii < *kv; ii++)
                AB[kk + ii] = 0.0;
        }
        AB[kk + *kv] = -1.0;
        AB[kk + *kv + 1] = 2.0;
        AB[kk + *kv + 2] = -1.0;
    }

    AB[0] = 0.0;
    if (*kv == 1)
        AB[1] = 0.0;
  
    AB[(*lab) * (*la) - 1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab,
                                           int *la, int *kv)
{
    for (size_t jj = 0; jj < (*la); jj++) {
        size_t kk = jj * (*lab);
        if (*kv >= 0) {
            for (size_t ii = 0; ii < *kv; ii++)
            	AB[kk + ii] = 0.0;
        }
        AB[kk + *kv]     = 0.0;
        AB[kk + *kv + 1] = 1.0;
        AB[kk + *kv + 2] = 0.0;
    }

    AB[1] = 0.0;
    AB[(*lab) * (*la) - 1] = 0.0;
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
    RHS[0] = *BC0;
    RHS[(*la) - 1] = *BC1;
    for (size_t jj = 1; jj < (*la) - 1; jj++)
        RHS[jj] = 0.0;
}  

void set_analytical_solution_DBC_1D(double *EX_SOL, double *X,
                                    int *la, double *BC0, double *BC1)
{
    double h; 
    double DELTA_T = (*BC1) - (*BC0);
    for (size_t jj = 0; jj < (*la); jj++)
        EX_SOL[jj] = (*BC0) + X[jj] * DELTA_T;
}  

void set_grid_points_1D(double *x, int *la)
{
    double h = 1.0 / (1.0 * ((*la) + 1));
    for (size_t jj = 0; jj < (*la); jj++)
        x[jj] = (jj + 1) * h;
}

void write_GB_operator_rowMajor_poisson1D(double *AB, int *lab,
                                          int *la, char *filename)
{
    FILE *file = fopen(filename, "wb");

    // Numbering from 1 to la
    if (file != NULL) {
        for (size_t ii = 0; ii < (*lab); ii++) {
            for (size_t jj = 0; jj < (*la); jj++)
                fprintf(file, "%lf\t", AB[ii * (*la) + jj]);
            fprintf(file, "\n");
        }
        fclose(file);
    } else {
        perror(filename);
    }
}

// TODO:
void write_GB_operator_colMajor_poisson1D(double *AB, int *lab,
                                          int *la, char *filename)
{
}

void write_vec(double *vec, int *la, char *filename)
{
    FILE *file = fopen(filename, "wb");

    // Numbering from 1 to la
    if (file != NULL) {
        for (size_t jj = 0; jj < (*la); jj++)
            fprintf(file, "%lf\n", vec[jj]);
        fclose(file);
    } else {
        perror(filename);
    } 
}  

void write_xy(double *vec, double *x, int *la, char *filename)
{
    FILE *file = fopen(filename, "wb");

    // Numbering from 1 to la
    if (file != NULL) {
        for (size_t jj = 0; jj < (*la); jj++)
            fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
        fclose(file);
    } else {
        perror(filename);
    } 
}  

void eig_poisson1D(double *eigval, int *la)
{
    for (size_t ii = 0; ii < *la; ii++) {
        double scal = (1.0 * ii + 1.0) * M_PI_2 * (1.0 / (*la + 1));
        eigval[ii] = 4 * sin(scal) * sin(scal);
    } 
}

double eigmax_poisson1D(int *la)
{
    double eigmax = sin(*la * M_PI_2 * (1.0 / (*la + 1)));
    return 4 * eigmax * eigmax;
}

double eigmin_poisson1D(int *la)
{
    double eigmin = sin(M_PI_2 * (1.0 / (*la + 1)));
    return 4 * eigmin * eigmin;
}

// TODO:
double richardson_alpha_opt(int *la)
{
}

// TODO:
void richardson_alpha(double *AB, double *RHS, double *X,
                      double *alpha_rich, int *lab, int *la,
                      int *ku, int*kl, double *tol, int *maxit)
{
}
