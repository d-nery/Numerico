/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "error.h"
#include "log.h"

matrix_t* matrix_create(int m, int n) {
    matrix_t* matrix;

    if (m <= 0 || n <= 0)
        error(ERR_NEGATIVE, "matrix_create");

    if ((matrix = (matrix_t*)calloc((size_t)1, sizeof(matrix_t))) == MAT_NULL)
        error(ERR_MEMORY, "matrix_create");

    matrix->l = m;
    matrix->c = n;

    if ((matrix->data = (double**)calloc(matrix->l, sizeof(double *))) == (double**)NULL) {
        free(matrix);
        error(ERR_MEMORY, "matrix_create");
    }

    for (int i = 0; i < matrix->l; i++)
        if ((matrix->data[i] = (double*)calloc((size_t)matrix->c, sizeof(double))) == (double*)NULL)
            error(ERR_MEMORY, "matrix_create");

    return matrix;
}

matrix_t* matrix_create_identity(int size) {
    matrix_t* matrix = matrix_create(size, size);

    for (int i = 0; i < size; i++)
        matrix_set(matrix, i, i, 1.0);

    return matrix;
}

matrix_t* matrix_copy(matrix_t* src, matrix_t* dst) {
    if (src == MAT_NULL)
        error(ERR_NULL, "matrix_copy");

    if (dst != MAT_NULL)
        matrix_free(dst);
    dst = matrix_create(src->l, src->c);

    for (int i = 0; i < src->l; i++)
        for (int j = 0; j < src->c; j++)
            matrix_set(dst, i, j, matrix_get(src, i, j));

    return dst;
}

void matrix_set(matrix_t* matrix, int i, int j, double v) {
    if (matrix == MAT_NULL)
        error(ERR_NULL, "matrix_set");

    if (i < 0 || i >= matrix->l || j < 0 || j >= matrix->c)
        error(ERR_OOB, "matrix_set");

    matrix->data[i][j] = v;
}

double matrix_get(const matrix_t* matrix, int i, int j) {
    if (matrix == MAT_NULL)
        error(ERR_NULL, "matrix_get");

    if (i < 0 || i >= matrix->l || j < 0 || j >= matrix->c)
        error(ERR_OOB, "matrix_get");

    return matrix->data[i][j];
}

void matrix_swap_lines(matrix_t* matrix, int l1, int l2) {
    if (matrix == MAT_NULL)
        error(ERR_NULL, "matrix_swap_lines");

    if (l1 < 0 || l1 >= matrix->l || l2 < 0 || l2 >= matrix->l)
        error(ERR_OOB, "matrix_swap_lines");

    double* aux = matrix->data[l1];
    matrix->data[l1] = matrix->data[l2];
    matrix->data[l2] = aux;
}

void print_matrix(matrix_t* matrix) {
    if (matrix == MAT_NULL)
        error(ERR_NULL, "print_matrix");

    if (log_get_level() <= LOG_INFO) {
        log_info("Matrix:");
        for (int i = 0; i < matrix->l; i++) {
            printf("|  ");
            for (int j = 0; j < matrix->c; j++) {
                printf("%+.3e  ", matrix_get(matrix, i, j));
            }
            printf("|\n");
        }
        printf("\n");
    }
}

matrix_t* matrix_add(matrix_t* A, matrix_t* B, matrix_t* R) {
    if (A == MAT_NULL || B == MAT_NULL)
        error(ERR_NULL, "matrix_add");

    if (A->l != B->l || A->c != B->c)
        error(ERR_SIZE, "matrix_add");

    if (R != MAT_NULL)
        matrix_free(R);

    R = matrix_create(A->l, A->c);

    for (int i = 0; i < A->l; i++)
        for (int j = 0; j < A->c; j++)
            matrix_set(R, i, j, matrix_get(A, i, j) + matrix_get(B, i, j));

    return R;
}

matrix_t* matrix_subtract(matrix_t* A, matrix_t* B, matrix_t* R) {
    return matrix_add(A, matrix_mult_scalar(-1.0, B, R), R);
}

matrix_t* matrix_mult_scalar(double n, matrix_t* A, matrix_t* R) {
    if (A == MAT_NULL)
        error(ERR_NULL, "matrix_set");

    if (R != MAT_NULL)
        matrix_free(R);

    R = matrix_create(A->l, A->c);

    for (int i = 0; i < A->l; i++)
        for (int j = 0; j < A->c; j++)
            matrix_set(R, i, j, matrix_get(A, i, j)*n);

    return R;
}

matrix_t* matrix_multiply(matrix_t* A, matrix_t* B, matrix_t* R) {
    if (A == MAT_NULL || B == MAT_NULL)
        error(ERR_NULL, "matrix_multiply");

    if (A->c != B->l)
        error(ERR_SIZE, "matrix_multiply");

    if (R != MAT_NULL)
        matrix_free(R);

    R = matrix_create(A->l, B->c);

    double sum = 0;
    for (int i = 0; i < R->l; i++) {
        for (int j = 0; j < R->c; j++) {
            sum = 0;
            for (int k = 0; k < A->c; k++) {
                sum += matrix_get(A, i, k) * matrix_get(B, k, j);
            }
            matrix_set(R, i, j, sum);
        }
    }

    return R;
}

void matrix_free(matrix_t* matrix) {
    if (matrix == MAT_NULL)
        error(ERR_NULL, "matrix_free");

    if (matrix->data == (double**)NULL)
        error(ERR_NULL, "matrix_free");

    for (int i = 0; i < matrix->l; i++) {
        if (matrix->data[i] == (double*)NULL)
            error(ERR_NULL, "matrix_free");
        free(matrix->data[i]);
    }

    matrix->l = 0;
    matrix->c = 0;
    free(matrix->data);
    free(matrix);
}
