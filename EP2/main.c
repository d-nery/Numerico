/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rkf.h"
#include "vector.h"
#include "matrix.h"

#define M_PI 3.141592653589793

#define min(a, b)          ((a) < (b) ? (a) : (b))
#define max(a, b)          ((a) > (b) ? (a) : (b))
#define constrain(v, a, b) ((v) < (a) ? (a) : (v) > (b) ? (b) : (v))

#define m_caso3 7

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

    matrix_t* A = matrix_create(m_caso3, m_caso3);
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

void calculo(vector_t* X0, double t0, double tf, double eps, double h, vector_t* f(double, vector_t*, vector_t*), char* name) {
    FILE* out = fopen(name, "w");

    vector_t** k = (vector_t**)calloc((size_t)6, sizeof(vector_t*));
    for (int i = 0; i < 6; i++)
        k[i] = vector_create(X0->size);

    double alpha;
    vector_t* xi   = vector_create(X0->size);
    vector_t* xi_b = vector_create(X0->size);

    double t = t0;
    vector_t* X = vector_copy(X0, X);
    vector_t* f_temp = VEC_NULL;
    vector_t* tal    = VEC_NULL;
    double maxtal = 999;
    int last = 0;

    while (t <= tf && last < 2) {
        maxtal = 999;
        while (maxtal > eps) {
            f_temp = f(t, X, f_temp);
            k[0] = vector_mult_scalar(h, f_temp, k[0]);

            f_temp = f(t + h_coeff[0]*h,
                vector_add(
                    X,
                    vector_mult_scalar(k2_coeff[0], k[0], NULL), NULL), f_temp);
            k[1] = vector_mult_scalar(h, f_temp, k[1]);

            f_temp = f(t + h_coeff[1]*h,
                vector_add_3(
                    X,
                    vector_mult_scalar(k3_coeff[0], k[0], NULL),
                    vector_mult_scalar(k3_coeff[1], k[1], NULL), NULL), f_temp);
            k[2] = vector_mult_scalar(h, f_temp, k[2]);

            f_temp = f(t + h_coeff[2]*h,
                vector_add_4(
                    X,
                    vector_mult_scalar(k4_coeff[0], k[0], NULL),
                    vector_mult_scalar(k4_coeff[1], k[1], NULL),
                    vector_mult_scalar(k4_coeff[2], k[2], NULL), NULL), f_temp);
            k[3] = vector_mult_scalar(h, f_temp, k[3]);

            f_temp = f(t + h_coeff[3]*h,
                vector_add_5(
                    X,
                    vector_mult_scalar(k5_coeff[0], k[0], NULL),
                    vector_mult_scalar(k5_coeff[1], k[1], NULL),
                    vector_mult_scalar(k5_coeff[2], k[2], NULL),
                    vector_mult_scalar(k5_coeff[3], k[3], NULL), NULL), f_temp);
            k[4] = vector_mult_scalar(h, f_temp, k[4]);

            f_temp = f(t + h_coeff[4]*h,
                vector_add_6(
                    X,
                    vector_mult_scalar(k6_coeff[0], k[0], NULL),
                    vector_mult_scalar(k6_coeff[1], k[1], NULL),
                    vector_mult_scalar(k6_coeff[2], k[2], NULL),
                    vector_mult_scalar(k6_coeff[3], k[3], NULL),
                    vector_mult_scalar(k6_coeff[4], k[4], NULL), NULL), f_temp);
            k[5] = vector_mult_scalar(h, f_temp, k[5]);

            xi   = vector_add_5(
                X,
                vector_mult_scalar(xi_coeff[0], k[0], NULL),
                vector_mult_scalar(xi_coeff[1], k[2], NULL),
                vector_mult_scalar(xi_coeff[2], k[3], NULL),
                vector_mult_scalar(xi_coeff[3], k[4], NULL), xi);

            xi_b = vector_add_6(
                X,
                vector_mult_scalar(xib_coeff[0], k[0], NULL),
                vector_mult_scalar(xib_coeff[1], k[2], NULL),
                vector_mult_scalar(xib_coeff[2], k[3], NULL),
                vector_mult_scalar(xib_coeff[3], k[4], NULL),
                vector_mult_scalar(xib_coeff[4], k[5], NULL), xi_b);
            tal = vector_abs(vector_subtract(xi_b, xi, tal), tal);
            tal = vector_mult_scalar(1.0/h, tal, tal);
            maxtal = vector_get(tal, 0);
            for (int i = 1; i < tal->size; i++)
                maxtal = max(maxtal, vector_get(tal, i));

            if (maxtal <= eps)
                break;

            alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
            h = alpha * h;
            h = min(h, tf - t);
            h = constrain(h, hmin, hmax);
        }
        t = t + h;
        if (t >= tf)
            last++;

        alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
        h = alpha * h;
        h = min(h, tf - t);
        h = constrain(h, hmin, hmax);
        X = vector_copy(xi, X);

        if (last < 2) {
            fprintf(out, "%.15e", t);
            for (int i = 0; i < xi->size; i++)
                fprintf(out, "  %.15e", vector_get(xi, i));
            fprintf(out, "\n");
        }
    }
    fclose(out);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s [caso]\n", argv[0]);
        exit(0);
    }

    vector_t* X0 = NULL;
    double t0    = 0.0;
    double tf    = 0.0;
    double eps   = 0.0;
    double h     = 0.0;

    switch(argv[1][0]) {
        case '1':
            X0 = vector_create(1);
            vector_set(X0, 0,  -18.95);
            t0  = 1.05;
            tf  = 3.0;
            eps = 1e-5;
            h   = 0.1;

            calculo(X0, t0, tf, eps, h, F1, "out.txt");
            break;

        case '2':
            X0 = vector_create(4);
            vector_set(X0, 0,  1);
            vector_set(X0, 1,  1);
            vector_set(X0, 2,  1);
            vector_set(X0, 3, -1);

            t0  = 0.0;
            tf  = 2.0;
            eps = 1e-5;
            h   = 0.1;

            calculo(X0, t0, tf, eps, h, F2, "out2.txt");
            break;

        case '3':
            X0 = vector_create(m_caso3);
            for (int i = 0; i < X0->size; i++)
                vector_set(X0, i, sin(M_PI*((double)(i+1)/((double)m_caso3 + 1.0))) + sin(m_caso3*M_PI*((double)(i + 1)/((double)m_caso3 + 1.0))));

            t0  = 0.0;
            tf  = 2.0;
            eps = 1e-5;
            h   = 0.1;

            calculo(X0, t0, tf, eps, h, F3, "out3.txt");
            break;

        default:
            printf("Erro! Escolha um caso valido\n");
    }
    return 0;
}
