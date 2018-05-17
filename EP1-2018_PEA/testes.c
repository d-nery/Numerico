/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#include <math.h>
#include <stdlib.h>

#define M_EULER 2.71828

#include "testes.h"
#include "error.h"

/**
 * Testes Iniciais propostos ao final do enunciado
 * FX e JFX são passadas para o metodo de Newton
 */

vector_t* F1(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "F1");

    vector_t* r = vector_create(x->size);

    double x1 = vector_get(x, 0);
    double x2 = vector_get(x, 1);

    vector_set(r, 0, 2*(x1 - 2));
    vector_set(r, 1, 3*(x2 - 3));

    return r;
}

matrix_t* JF1(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "JF1");

    matrix_t* J = matrix_create(x->size, x->size);

    double x1 = vector_get(x, 0);
    double x2 = vector_get(x, 1);

    matrix_set(J, 0, 0, 2);
    matrix_set(J, 0, 1, 0);

    matrix_set(J, 1, 0, 0);
    matrix_set(J, 1, 1, 3);

    return J;
}


vector_t* F2(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "F2");

    vector_t* r = vector_create(x->size);

    double x1 = vector_get(x, 0);
    double x2 = vector_get(x, 1);
    double x3 = vector_get(x, 2);
    double x4 = vector_get(x, 3);

    vector_set(r, 0,  4*x1 - 1*x2 + 1*x3 - x1*x4);
    vector_set(r, 1, -1*x1 + 3*x2 - 2*x3 - x2*x4);
    vector_set(r, 2,  1*x1 - 2*x2 + 3*x3 - x3*x4);
    vector_set(r, 3,  x1*x1 + x2*x2 + x3*x3 - 1);

    return r;
}

matrix_t* JF2(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "JF2");

    matrix_t* J = matrix_create(x->size, x->size);

    double x1 = vector_get(x, 0);
    double x2 = vector_get(x, 1);
    double x3 = vector_get(x, 2);
    double x4 = vector_get(x, 3);

    matrix_set(J, 0, 0,  4 - x4);
    matrix_set(J, 0, 1, -1);
    matrix_set(J, 0, 2,  1);
    matrix_set(J, 0, 3, -x1);

    matrix_set(J, 1, 0, -1);
    matrix_set(J, 1, 1, 3 - x4);
    matrix_set(J, 1, 2, -2);
    matrix_set(J, 1, 3, -x2);

    matrix_set(J, 2, 0, 1);
    matrix_set(J, 2, 1, -2);
    matrix_set(J, 2, 2, 3 - x4);
    matrix_set(J, 2, 3, -x3);

    matrix_set(J, 3, 0, 2*x1);
    matrix_set(J, 3, 1, 2*x2);
    matrix_set(J, 3, 2, 2*x3);
    matrix_set(J, 3, 3, 0);

    return J;
}

vector_t* F3(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "F3");

    const int n = x->size + 1;

    vector_t* r = vector_create(n - 1);

    vector_set(r, 0, 2*vector_get(x, 0) - vector_get(x, 1) - pow(M_EULER, vector_get(x, 0))/pow(n, 2));
    for (int i = 1; i < n-2; i++) {
        vector_set(r, i, -vector_get(x, i-1) + 2*vector_get(x, i) - vector_get(x, i + 1) - pow(M_EULER, vector_get(x, i))/pow(n, 2));
    }
    vector_set(r, n-2, -vector_get(x, n-3) + 2*vector_get(x, n-2) - pow(M_EULER, vector_get(x, n-2))/pow(n, 2));

    return r;
}

matrix_t* JF3(vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "JF3");

    const int n = x->size + 1;

    matrix_t* J = matrix_create(n - 1, n - 1);

    matrix_set(J, 0, 0, 2 - pow(M_EULER, vector_get(x, 0))/pow(n, 2));
    matrix_set(J, 0, 1, -1);

    for (int i = 1; i < n - 2; i++) {
            matrix_set(J, i, i - 1, -1);
            matrix_set(J, i, i, 2 - pow(M_EULER, vector_get(x, i))/pow(n, 2));
            matrix_set(J, i, i + 1, -1);
    }

    matrix_set(J, n - 2, n - 3, -1);
    matrix_set(J, n - 2, n - 2, 2 - pow(M_EULER, vector_get(x, n - 2))/pow(n, 2));

    return J;
}
