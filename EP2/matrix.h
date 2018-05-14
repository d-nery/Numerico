/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
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

void print_matrix(matrix_t* matrix);

void matrix_set(matrix_t* matrix, int i, int j, double v);
double matrix_get(const matrix_t* matrix, int i, int j);

// Operations
matrix_t* matrix_add(matrix_t* A, matrix_t* B);
matrix_t* matrix_subtract(matrix_t* A, matrix_t* B);
matrix_t* matrix_mult_scalar(double n, matrix_t* A);

void matrix_free(matrix_t* matrix);

#endif
