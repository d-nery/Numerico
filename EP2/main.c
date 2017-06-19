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

#ifndef CASE
#error "define CASE in makeflie"
#endif

#define min(a, b)          ((a) < (b) ? (a) : (b))
#define max(a, b)          ((a) > (b) ? (a) : (b))
#define constrain(v, a, b) ((v) < (a) ? (a) : (v) > (b) ? (b) : (v))

double F(double t, double x) {
    return 1.0 + pow(x - t, 2);
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

// double F2(double t, matrix_t* X) {

//     return matrix_multiply(A, X);
// }

void calculo1(double x0, double t0, double tf, double eps, double h, double f(double, double), char* name) {
    FILE* out = fopen(name, "w");

    double xi, xi_b, alpha;
    double tal = 999;

    double k[6] = { 0 };

    double t = t0;
    double x = x0;

    int last = 0;

    while (t <= tf && last < 2) {
        printf("-- Iteracao t: %f\n", t);
        tal = 999;
        while (tal > eps) {
            k[0] = h * f(t, x);
            k[1] = h * f(t + h_coeff[0]*h, x + k2_coeff[0]*k[0]);
            k[2] = h * f(t + h_coeff[1]*h, x + k3_coeff[0]*k[0] + k3_coeff[1]*k[1]);
            k[3] = h * f(t + h_coeff[2]*h, x + k4_coeff[0]*k[0] + k4_coeff[1]*k[1] + k4_coeff[2]*k[2]);
            k[4] = h * f(t + h_coeff[3]*h, x + k5_coeff[0]*k[0] + k5_coeff[1]*k[1] + k5_coeff[2]*k[2] + k5_coeff[3]*k[3]);
            k[5] = h * f(t + h_coeff[4]*h, x + k6_coeff[0]*k[0] + k6_coeff[1]*k[1] + k6_coeff[2]*k[2] + k6_coeff[3]*k[3] + k6_coeff[4]*k[4]);

            xi   = x + xi_coeff[0]*k[0]  + xi_coeff[1]*k[2]  + xi_coeff[2]*k[3]  + xi_coeff[3]*k[4];
            xi_b = x + xib_coeff[0]*k[0] + xib_coeff[1]*k[2] + xib_coeff[2]*k[3] + xib_coeff[3]*k[4] + xib_coeff[4]*k[5];

            tal  = fabs(xi_b - xi)/h;
            if (tal <= eps)
                break;

            alpha = pow((eps)/(c_security*tal), 1.0/4.0);
            h = alpha * h;
            h = min(h, tf - t);
            h = constrain(h, hmin, hmax);
        }
        t = t + h;
        if (t >= tf)
            last++;

        alpha = pow((eps)/(c_security*fabs(tal)), 1.0/4.0);
        h = alpha * h;
        h = min(h, tf - t);
        h = constrain(h, hmin, hmax);
        x = xi;

        printf("%.6f\n%.6f\n%.6f\n%.6f\n%.6f\n%.6f\n%.6f\n%.6f\n%.6e\n%.5f\n\n", k[0], k[1], k[2], k[3], k[4], k[5], xi, xi_b, tal, h);
        if (last < 2)
            fprintf(out, "%.6f  %.6f\n", t, xi);
    }
    fclose(out);
}

void calculo2(vector_t* X0, double t0, double tf, double eps, double h, vector_t* f(double, vector_t*, vector_t*), char* name) {
    vector_t** k = (vector_t**)calloc((size_t)6, sizeof(vector_t*));
    for (int i = 0; i < 6; i++)
        k[i] = vector_create(X0->size);

    double alpha;
    vector_t* xi   = vector_create(X0->size);
    vector_t* xi_b = vector_create(X0->size);

    double t = t0;
    vector_t* X = X0;
    vector_t* f_temp = vector_create(X0->size);
    vector_t* tal;
    double maxtal = 999;

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
        printf("Maxtal: %f\n", maxtal);

        alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
        h = alpha * h;
        h = min(h, tf - t);
        h = constrain(h, hmin, hmax);
        printf("Prox h: %f", h);

        printf("\nk0:");
        print_vector(k[0]);
        printf("\nk1:");
        print_vector(k[1]);
        printf("\nk2:");
        print_vector(k[2]);
        printf("\nk3:");
        print_vector(k[3]);
        printf("\nk4:");
        print_vector(k[4]);
        printf("\nk5:");
        print_vector(k[5]);
        printf("\nxi:");
        print_vector(xi);
        printf("\nxib:");
        print_vector(xi_b);
        printf("\ntal:");
        print_vector(tal);
        printf("\n");
        printf("\n");
    }
}

int main() {
#if CASE == 1
    double x0  = -18.95;
    double t0  = 1.05;
    double tf  = 3.0;
    double eps = 1e-5;
    double h   = 0.1;

    calculo1(x0, t0, tf, eps, h, F, "out.txt");
#elif CASE == 2
    vector_t* X0 = vector_create(4);
    vector_set(X0, 0,  1);
    vector_set(X0, 1,  1);
    vector_set(X0, 2,  1);
    vector_set(X0, 3, -1);

    double t0  = 0.0;
    double tf  = 2.0;
    double eps = 1e-5;
    double h   = 0.1;

    calculo2(X0, t0, tf, eps, h, F2, "out2.txt");
#else
#error CASE nao encontrado
#endif
    return 0;
}
