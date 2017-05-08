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

	for (int i = 0; i < matrix->l; i++) {
		if ((matrix->data[i] = (double*)calloc((size_t)matrix->c, sizeof(double))) == (double*)NULL)
			error(ERR_MEMORY, "matrix_create");
	}

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

		matrix->data[l][c] = val;
	}

	fclose(fp);
	return matrix;
}

vector_t* vector_create(int size) {
	return NULL;
}

void print_matrix(matrix_t* matrix) {
	for (int i = 0; i < matrix->l; i++) {
		for (int j = 0; j < matrix->c; j++) {
			if (matrix->data[i][j] >= 0)
				printf(" ");
			printf("%.12e ", matrix->data[i][j]);
		}
		printf("\n");
	}
}
