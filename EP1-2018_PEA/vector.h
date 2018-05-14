/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include "matrix.h"

#define VEC_NULL (vector_t *)NULL

#define sgn(v) ((v) < 0 ? -1 : 1)

typedef struct {
    double* data;
    int size;
} vector_t;

vector_t* vector_create(const int size);
vector_t* vector_create_from_file(const char* name, const int lines);
vector_t* vector_copy(const vector_t* v, vector_t* r);

void print_vector(const vector_t* vector);
void output_vector(const vector_t* vector, const char* filename);

void vector_set(vector_t* vector, const int i, const double v);
double vector_get(const vector_t* vector, const int i);

double vector_norm(const vector_t* vector);
double vector_multiply(const vector_t* u, const vector_t* v);

vector_t* vector_add(const vector_t* u, const vector_t* v, vector_t* r);

vector_t* vector_subtract(const vector_t* u, const vector_t* v, vector_t* r);
vector_t* vector_mult_scalar(const double n, const vector_t* u, vector_t* r);
vector_t* vector_mult_matrix(const matrix_t* A, const vector_t* v, vector_t* r);

vector_t* vector_abs(const vector_t* vector, vector_t* r);

void vector_free(vector_t* vector);

#endif
