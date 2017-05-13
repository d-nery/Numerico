/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#define VEC_NULL (vector_t *)NULL

typedef struct {
	double* data;
	int size;
} vector_t;

vector_t* vector_create(int size);
vector_t* vector_create_from_file(char* name, int lines);
void print_vector(vector_t* vector);

void vector_set(vector_t* vector, int i, float v);
float vector_get(vector_t* vector, int i);

vector_t* vector_add(vector_t* A, vector_t* B);
vector_t* vector_subtract(vector_t* A, vector_t* B);
vector_t* vector_mult_scalar(float n, vector_t* A);

#endif
