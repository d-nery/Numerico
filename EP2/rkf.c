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

#include "rkf.h"

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
    int last = 0;

    vector_t* xi   = vector_create(X0->size);
    vector_t* xi_b = vector_create(X0->size);

    vector_t* X = vector_copy(X0, X);
    vector_t* f_temp = VEC_NULL;
    vector_t* tal    = VEC_NULL;

    vector_t* mult_temp = VEC_NULL;
    vector_t* add_temp  = VEC_NULL;

    int it = 0;
    int go = 0;
    while (t <= tf && last < 2 && !go) {
        printf("Iteracao: %d\n", it++);
        printf("-- t: %.12lf\n", t);
        printf("-- h: %.12lf\n", h);
        maxtal = 999;
        while (maxtal > eps) {
            f_temp = f(t, X, f_temp);
            k[0] = vector_mult_scalar(h, f_temp, k[0]);

            mult_temp = vector_mult_scalar(k2_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            f_temp = f(t + h_coeff[0]*h, add_temp, f_temp);
            k[1] = vector_mult_scalar(h, f_temp, k[1]);

            mult_temp = vector_mult_scalar(k3_coeff[0], k[0], mult_temp);
            add_temp  = vector_add(X, mult_temp, add_temp);
            mult_temp = vector_mult_scalar(k3_coeff[1], k[1], mult_temp);
            add_temp  = vector_add(add_temp, mult_temp, add_temp);
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

            // for (int i = 0; i < 6; i++) {
            //     printf("\n-- k%d:\n", i+1);
            //     print_vector(k[i]);
            // }
            // printf("-- xi:\n");
            // print_vector(xi);
            // printf("-- xb:\n");
            // print_vector(xi_b);
            // printf("-- tal:\n");
            // print_vector(tal);
            // printf("-- maxtal: %.12lf\n", maxtal);
            // printf("-- alpha:  %.12lf\n", alpha);
            // printf("-- prox h: %.12lf\n", h);
            // go = 1;
            // break;
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

void rkf45_error(double xt(double), char* file) {

}
