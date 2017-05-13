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

#include "vector.h"
#include "error.h"

void print_vector(vector_t* vector) {
	if (vector == VEC_NULL)
		error(ERR_NULL, "print_vector");

	for (int i = 0; i < vector->size; i++) {
		if (vector->data[i] >= 0)
			printf(" ");
		printf("%.12e\n", vector->data[i]);
	}
}

vector_t* vector_create(int size) {
	vector_t* vector;

	if (size <= 0)
		error(ERR_NEGATIVE, "vector_create");

	if ((vector = (vector_t*)calloc((size_t)1, sizeof(vector_t))) == VEC_NULL)
		error(ERR_MEMORY, "vector_create");

	vector->size = size;

	if ((vector->data = (double*)calloc(vector->size, sizeof(double))) == (double*)NULL) {
		free(vector);
		error(ERR_MEMORY, "vector_create");
	}

	return vector;
}

vector_t* vector_create_from_file(char* name, int lines) {
	vector_t* vector;
	FILE* fp = fopen(name, "r");

	if (fp == NULL)
		error(ERR_OPEN_FILE, "vector_create_from_file");

	vector = vector_create(lines);

	double val;
	for (int i = 0; i < lines; i++) {
		if (fscanf(fp, "%lf", &val) != 1)
			error(ERR_READ_FILE, "vector_create_from_file");

		vector->data[i] = val;
	}

	fclose(fp);
	return vector;
}

void vector_set(vector_t* vector, int i, float v) {
	if (vector == VEC_NULL)
		error(ERR_NULL, "vector_set");

	if (i < 0 || i > vector->size)
		error(ERR_OOB, "vector_set");

	vector->data[i] = v;
}

float vector_get(vector_t* vector, int i) {
	if (vector == VEC_NULL)
		error(ERR_NULL, "vector_get");

	if (i < 0 || i > vector->size)
		error(ERR_OOB, "vector_get");

	return vector->data[i];
}


vector_t* vector_add(vector_t* u, vector_t* v) {
	if (u == VEC_NULL || v == VEC_NULL)
		error(ERR_NULL, "vector_set");

	if (u->size != v->size)
		error(ERR_SIZE, "vector_add");

	vector_t* result;
	result = vector_create(u->size);

	for (int i = 0; i < u->size; i++)
		vector_set(result, i, vector_get(u, i) + vector_get(v, i));

	return result;
}
