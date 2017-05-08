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
#define VEC_NULL (vector_t *)NULL

typedef struct {
	double** data;
	int l, c;
} matrix_t;

typedef struct {
	double* data;
	int size;
} vector_t;

matrix_t* matrix_create(int m, int n);
matrix_t* matrix_create_from_file(char* name);
vector_t* vector_create(int size);
vector_t* vector_create_from_file(char* name);

void print_matrix(matrix_t* matrix);

#endif
