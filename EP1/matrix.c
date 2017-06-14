/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"
#include "error.h"

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

matrix_t* matrix_create_from_file(char* name) {
	matrix_t* matrix;
	FILE* fp = fopen(name, "r");
	int m = 0, n = 0, lines;

	if (fp == NULL)
		error(ERR_OPEN_FILE, "matrix_create_from_file");

	if (fscanf(fp, "%d %d %d", &m, &n, &lines) < 3)
		error(ERR_READ_FILE, "matrix_create_from_file");

	matrix = matrix_create(m, n);

	int l, c;
	double val;
	for (int i = 0; i < lines; i++) {
		if (fscanf(fp, "%d %d %lf", &l, &c, &val) < 3)
			error(ERR_READ_FILE, "matrix_create_from_file");

		matrix_set(matrix, l, c, val);
	}

	fclose(fp);
	return matrix;
}

void matrix_set(matrix_t* matrix, int i, int j, double v) {
	if (matrix == MAT_NULL)
		error(ERR_NULL, "matrix_set");

	if (i < 0 || i > matrix->l || j < 0 || j > matrix->c)
		error(ERR_OOB, "matrix_set");

	matrix->data[i][j] = v;
}

double matrix_get(matrix_t* matrix, int i, int j) {
	if (matrix == MAT_NULL)
		error(ERR_NULL, "matrix_get");

	if (i < 0 || i > matrix->l || j < 0 || j > matrix->c)
		error(ERR_OOB, "matrix_get");

	return matrix->data[i][j];
}

void print_matrix(matrix_t* matrix) {
	if (matrix == MAT_NULL)
		error(ERR_NULL, "print_matrix");

	for (int i = 0; i < matrix->l; i++) {
		for (int j = 0; j < matrix->c; j++) {
			if (matrix_get(matrix, i, j) >= 0)
				printf(" ");
			printf("%.3e ", matrix_get(matrix, i, j));
		}
		printf("\n");
	}
}

matrix_t* matrix_add(matrix_t* A, matrix_t* B) {
	if (A == MAT_NULL || B == MAT_NULL)
		error(ERR_NULL, "matrix_set");

	if (A->l != B->l || A->c != B->c)
		error(ERR_SIZE, "matrix_add");

	matrix_t* result;
	result = matrix_create(A->l, A->c);

	for (int i = 0; i < A->l; i++)
		for (int j = 0; j < A->c; j++)
			matrix_set(result, i, j, matrix_get(A, i, j) + matrix_get(B, i, j));

	return result;
}

matrix_t* matrix_subtract(matrix_t* A, matrix_t* B) {
	return matrix_add(A, matrix_mult_scalar(-1.0, B));
}

matrix_t* matrix_mult_scalar(double n, matrix_t* A) {
	if (A == MAT_NULL)
		error(ERR_NULL, "matrix_set");

	matrix_t* result;
	result = matrix_create(A->l, A->c);

	for (int i = 0; i < A->l; i++)
		for (int j = 0; j < A->c; j++)
			matrix_set(result, i, j, matrix_get(A, i, j)*n);

	return result;
}
