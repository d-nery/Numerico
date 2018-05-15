/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

#define MAT_NULL (matrix_t *)NULL

typedef struct {
    double** data;
    int l, c;
} matrix_t;

matrix_t* matrix_create(int m, int n);
matrix_t* matrix_create_identity(int size);
matrix_t* matrix_create_from_file(char* name);
matrix_t* matrix_copy(matrix_t* src, matrix_t* dst);

void print_matrix(matrix_t* matrix);

void matrix_set(matrix_t* matrix, int i, int j, double v);
double matrix_get(const matrix_t* matrix, int i, int j);

void matrix_swap_lines(matrix_t* matrix, int l1, int l2);

// Operations
matrix_t* matrix_add(matrix_t* A, matrix_t* B, matrix_t* R);
matrix_t* matrix_subtract(matrix_t* A, matrix_t* B, matrix_t* R);
matrix_t* matrix_mult_scalar(double n, matrix_t* A, matrix_t* R);
matrix_t* matrix_multiply(matrix_t* A, matrix_t* B, matrix_t* R);

void matrix_free(matrix_t* matrix);

#endif
