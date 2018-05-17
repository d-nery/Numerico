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
#include <math.h>

#include "vector.h"
#include "error.h"
#include "log.h"

void print_vector(const vector_t* vector) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "print_vector");

    if (log_get_level() <= LOG_INFO) {
        log_info("Vector:");
        for (int i = 0; i < vector->size; i++) {
            if (vector_get(vector, i) >= 0)
                printf(" ");
            printf("%.5e ", vector_get(vector, i));
            if (vector->size > 9)
                printf("\n");
        }
        printf("\n");
    }
}

vector_t* vector_create(const int size) {
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

vector_t* vector_create_from_file(const char* name, const int lines) {
    vector_t* vector = VEC_NULL;
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

vector_t* vector_copy(const vector_t* v, vector_t* r) {
    if (v == VEC_NULL)
        error(ERR_NULL, "vector_abs");

    if (r == VEC_NULL)
        r = vector_create(v->size);

    if (r->size != v->size) {
        vector_free(r);
        r = vector_create(v->size);
    }

    for (int i = 0; i < v->size; i++)
        vector_set(r, i, vector_get(v, i));

    return r;
}

void vector_set(vector_t* vector, const int i, const double v) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "vector_set");

    if (i < 0 || i > vector->size)
        error(ERR_OOB, "vector_set");

    vector->data[i] = v;
}

double vector_get(const vector_t* vector, const int i) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "vector_get");

    if (i < 0 || i > vector->size)
        error(ERR_OOB, "vector_get");

    return vector->data[i];
}

double vector_norm(const vector_t* vector) {
    return sqrt(vector_multiply(vector, vector));
}

double vector_multiply(const vector_t* u, const vector_t* v) {
    if (u == VEC_NULL || v == VEC_NULL)
        error(ERR_NULL, "vector_multiply");

    if (u->size != v->size)
        error(ERR_SIZE, "vector_multiply");

    double result = 0;

    for (int i = 0; i < u->size; i++)
        result += vector_get(u, i) * vector_get(v, i);

    return result;
}

vector_t* vector_add(const vector_t* u, const vector_t* v, vector_t* r) {
    if (u == VEC_NULL || v == VEC_NULL)
        error(ERR_NULL, "vector_set");

    if (u->size != v->size)
        error(ERR_SIZE, "vector_add");

    if (r == VEC_NULL)
        r = vector_create(u->size);

    if (r->size != u->size) {
        vector_free(r);
        r = vector_create(u->size);
    }

    for (int i = 0; i < u->size; i++)
        vector_set(r, i, vector_get(u, i) + vector_get(v, i));

    return r;
}

vector_t* vector_subtract(const vector_t* u, const vector_t* v, vector_t* r) {
    if (u == VEC_NULL || v == VEC_NULL)
        error(ERR_NULL, "vector_subtract");

    if (u->size != v->size)
        error(ERR_SIZE, "vector_subtract");

    if (r == VEC_NULL)
        r = vector_create(u->size);

    if (r->size != u->size) {
        vector_free(r);
        r = vector_create(u->size);
    }

    for (int i = 0; i < u->size; i++)
        vector_set(r, i, vector_get(u, i) - vector_get(v, i));

    return r;
}

vector_t* vector_mult_scalar(const double n, const vector_t* u, vector_t* r) {
    if (u == VEC_NULL)
        error(ERR_NULL, "vector_mult_scalar");

    if (r == VEC_NULL)
        r = vector_create(u->size);

    if (r->size != u->size) {
        vector_free(r);
        r = vector_create(u->size);
    }

    for (int i = 0; i < u->size; i++)
        vector_set(r, i, vector_get(u, i)*n);

    return r;
}

void vector_free(vector_t* vector) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "vector_free");

    if (vector->data == (double*)NULL)
        error(ERR_NULL, "vector_free");

    vector->size = 0;
    free(vector->data);
    free(vector);
}

void output_vector(const vector_t* vector, const char* filename) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "print_vector");

    FILE* f = fopen(filename, "w");

    for (int i = 0; i < vector->size; i++) {
        if (vector_get(vector, i) >= 0)
            fprintf(f, " ");
        fprintf(f,  "%.12e\n", vector_get(vector, i));
    }
}

vector_t* vector_mult_matrix(const matrix_t* A, const vector_t* v, vector_t* r) {
    if (A == MAT_NULL || v == VEC_NULL)
        error(ERR_NULL, "vector_mult_matrix");

    if (A->data == (double**)NULL || v->data == (double*)NULL)
        error(ERR_NULL, "vector_mult_matrix");

    if (A->c != v->size)
        error(ERR_SIZE, "vector_mult_matrix");

    vector_free(r);
    r = vector_create(A->l);

    for (int i = 0; i < A->l; i++) {
        for (int j = 0; j < A->c; j++) {
            vector_set(r, i, vector_get(r, i) + vector_get(v, j)*matrix_get(A, i, j));
        }
    }

    return r;
}

vector_t* vector_abs(const vector_t* vector, vector_t* r) {
    if (vector == VEC_NULL)
        error(ERR_NULL, "vector_abs");

    if (r == VEC_NULL)
        r = vector_create(vector->size);

    if (r->size != vector->size) {
        vector_free(r);
        r = vector_create(vector->size);
    }

    for (int i = 0; i < vector->size; i++)
        vector_set(r, i, fabs(vector_get(vector, i)));

    return r;
}
