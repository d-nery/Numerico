/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include "rkf.h"

// #define DEBUG

// Coeficientes do RKF45
static const double h_coeff[]   = { 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0 };
static const double k2_coeff[]  = { 1.0/4.0 };
static const double k3_coeff[]  = { 3.0/32.0, 9.0/32.0 };
static const double k4_coeff[]  = { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0 };
static const double k5_coeff[]  = { 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0 };
static const double k6_coeff[]  = { -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0 };
static const double xi_coeff[]  = { 25.0/216.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0 };
static const double xib_coeff[] = { 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 };

void rkf45_solve(vector_t* X0, double t0, double tf, double eps, double h, vector_t* f(double, vector_t*, vector_t*), char* name) {
    // Arquivo para os dados de saida
    FILE* out = fopen(name, "w");

    vector_t** k = (vector_t**)calloc((size_t)6, sizeof(vector_t*));
    for (int i = 0; i < 6; i++)
        k[i] = vector_create(X0->size);

    double alpha;

    double t = t0;
    double maxtal = 999;
    // int last = 0;

    vector_t* xi   = vector_create(X0->size);
    vector_t* xi_b = vector_create(X0->size);

    vector_t* X = vector_copy(X0, X);
    vector_t* f_temp = VEC_NULL;
    vector_t* tal    = VEC_NULL;

    vector_t* mult_temp = VEC_NULL;
    vector_t* add_temp  = VEC_NULL;

    int it = 0;
    while (t < tf) {
        printf("Iteracao: %d\n", it++);
        printf("-- t: %.12lf\n", t);
        printf("-- h: %.12lf\n", h);
        maxtal = eps + 1;
        while (maxtal > eps) {
            f_temp = f(t, X, f_temp);
            k[0] = vector_mult_scalar(h, f_temp, k[0]);

            mult_temp = vector_mult_scalar(k2_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            f_temp = f(t + h_coeff[0]*h, add_temp, f_temp);
            k[1] = vector_mult_scalar(h, f_temp, k[1]);

            mult_temp = vector_mult_scalar(k3_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 1; i++) {
                mult_temp = vector_mult_scalar(k3_coeff[i], k[i], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            f_temp = f(t + h_coeff[1]*h, add_temp, f_temp);
            k[2] = vector_mult_scalar(h, f_temp, k[2]);

            mult_temp = vector_mult_scalar(k4_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 2; i++) {
                mult_temp = vector_mult_scalar(k4_coeff[i], k[i], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            f_temp = f(t + h_coeff[2]*h, add_temp, f_temp);
            k[3] = vector_mult_scalar(h, f_temp, k[3]);

            mult_temp = vector_mult_scalar(k5_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 3; i++) {
                mult_temp = vector_mult_scalar(k5_coeff[i], k[i], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            f_temp = f(t + h_coeff[3]*h, add_temp, f_temp);
            k[4] = vector_mult_scalar(h, f_temp, k[4]);

            mult_temp = vector_mult_scalar(k6_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 4; i++) {
                mult_temp = vector_mult_scalar(k6_coeff[i], k[i], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            f_temp = f(t + h_coeff[4]*h, add_temp, f_temp);
            k[5] = vector_mult_scalar(h, f_temp, k[5]);

            mult_temp = vector_mult_scalar(xi_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 3; i++) {
                mult_temp = vector_mult_scalar(xi_coeff[i], k[i+1], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            xi   = vector_copy(add_temp, xi);

            mult_temp = vector_mult_scalar(xib_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            for (int i = 1; i <= 4; i++) {
                mult_temp = vector_mult_scalar(xib_coeff[i], k[i+1], mult_temp);
                add_temp  = vector_add(add_temp, mult_temp, add_temp);
            }
            xi_b = vector_copy(add_temp, xi_b);

            tal = vector_abs(vector_subtract(xi_b, xi, tal), tal);
            maxtal = vector_get(tal, 0);
            for (int i = 1; i < tal->size; i++)
                maxtal = max(maxtal, vector_get(tal, i));
            maxtal /= h;

#ifdef DEBUG
printf("k1: "); print_vector(k[0]);
printf("k2: "); print_vector(k[1]);
printf("k3: "); print_vector(k[2]);
printf("k4: "); print_vector(k[3]);
printf("k5: "); print_vector(k[4]);
printf("k6: "); print_vector(k[5]);
printf("x(i+1):  "); print_vector(xi);
printf("x(i+1)_: "); print_vector(xi_b);
printf("Tau maximo: %.8e\n", maxtal);
#endif
            if (maxtal <= eps) {
#ifdef DEBUG
printf("Resposta aceita\n");
#endif
                break;
            } else {
#ifdef DEBUG
printf("Resposta rejeitada\n");
#endif
            }

            alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
            h = alpha * h;
            h = constrain(h, hmin, hmax);
            h = min(h, tf - t);
#ifdef DEBUG
printf("Alfa: %.8f\n", alpha);
printf("Novo h: %.8f\n\n", h);
#endif

            usleep(100*1000);
        }
        t = t + h;
        X = vector_copy(xi, X);
// printf("ti: %.8f  h: %.8f\n", t, h);

        alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
// printf("Alfa: %.8f\n", alpha);
        h = alpha * h;
        h = constrain(h, hmin, hmax);
        h = min(h, tf - t);
// printf("Novo h: %.8f\n\n", h);

        fprintf(out, "%.15e", t);
        for (int i = 0; i < xi->size; i++)
            fprintf(out, "  %.15e", vector_get(xi, i));
        fprintf(out, "\n");
    }
    fclose(out);
}

vector_t* F1(double t, vector_t* X, vector_t* res) {
    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    for (int i = 0; i < X->size; i++)
        vector_set(res, i, 1.0 + pow(vector_get(X, i) - t, 2.0));

    return res;
}

vector_t* F2(double t, vector_t* X, vector_t* res) {
    (void)t; // Avoids unused warning

    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    matrix_t* A = matrix_create(4, 4);
    matrix_set(A, 0, 0, -2);
    matrix_set(A, 0, 1, -1);
    matrix_set(A, 0, 2, -1);
    matrix_set(A, 0, 3, -2);
    matrix_set(A, 1, 0,  1);
    matrix_set(A, 1, 1, -2);
    matrix_set(A, 1, 2,  2);
    matrix_set(A, 1, 3, -1);
    matrix_set(A, 2, 0, -1);
    matrix_set(A, 2, 1, -2);
    matrix_set(A, 2, 2, -2);
    matrix_set(A, 2, 3, -1);
    matrix_set(A, 3, 0,  2);
    matrix_set(A, 3, 1, -1);
    matrix_set(A, 3, 2,  1);
    matrix_set(A, 3, 3, -2);

    res = vector_mult_matrix(A, X, res);

    matrix_free(A);

    return res;
}

vector_t* F3(double t, vector_t* X, vector_t* res) {
    (void)t; // Avoids unused warning

    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    matrix_t* A = matrix_create(X->size, X->size);
    for (int i = 0; i < A->l; i++) {
        for (int j = 0; j < A->c; j++) {
            if (i == j)
                matrix_set(A, i, j, -2.0);
            else if ((i < A->l - 1 && j == i + 1) || (j < A->c - 1 && i == j + 1))
                matrix_set(A, i, j, 1);
        }
    }

    res = vector_mult_matrix(A, X, res);

    matrix_free(A);

    return res;
}

vector_t* F_chua(double t, vector_t* X, vector_t* res) {
    (void)t; // Avoids unused warning

    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    matrix_t* A = matrix_create(3, 3);
    matrix_set(A, 0, 0, -1.0/(R*C1));
    matrix_set(A, 0, 1,  1.0/(R*C1));
    matrix_set(A, 0, 2,  0.0);
    matrix_set(A, 1, 0,  1.0/(R*C2));
    matrix_set(A, 1, 1, -1.0/(R*C2));
    matrix_set(A, 1, 2, -1.0/(C2));
    matrix_set(A, 2, 0,  0.0);
    matrix_set(A, 2, 1, -1.0/(L));
    matrix_set(A, 2, 2,  0.0);

    res = vector_mult_matrix(A, X, res);

	vector_set(res, 0, vector_get(res, 0) - g(vector_get(X, 0))/C1);

    matrix_free(A);

    return res;
}
