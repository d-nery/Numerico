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
#include "chua.h"

#define DEBUG 1

// Coeficientes do RKF45
static const double h_coeff[]    = { 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0 };

static const double k_coeff[][5] = {{ 0.0,            0.0,           0.0,           0.0,            0.0       },
                                    { 1.0/4.0,        0.0,           0.0,           0.0,            0.0       },
                                    { 3.0/32.0,       9.0/32.0,      0.0,           0.0,            0.0       },
                                    { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0,            0.0       },
                                    { 439.0/216.0,   -8.0,           3680.0/513.0, -845.0/4104.0,   0.0       },
                                    { -8.0/27.0,      2.0,          -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0 }};

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
        maxtal = eps + 1;
        while (maxtal > eps) {
            // Calculo dos ks
            for (int j = 0; j < 6; j++) {
                mult_temp = vector_mult_scalar(k_coeff[j][0], k[0], mult_temp);
                add_temp  = vector_add(X, mult_temp, add_temp);
                for (int i = 1; i <= j-1; i++) {
                    mult_temp = vector_mult_scalar(k_coeff[j][i], k[i], mult_temp);
                    add_temp  = vector_add(add_temp, mult_temp, add_temp);
                }
                f_temp = f(t + h_coeff[j]*h, add_temp, f_temp);
                k[j] = vector_mult_scalar(h, f_temp, k[j]);
            }

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

#if DEBUG > 2
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
#if DEBUG > 1
                printf("Resposta aceita\n");
#endif
                break;
            } else {
#if DEBUG > 1
                printf("Resposta rejeitada\n");
#endif
            }

            alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
            h = alpha * h;
            // h = constrain(h, hmin, hmax);
            // h = min(h, tf - t);
#if DEBUG > 1
            printf("Alfa:   %.8e\n", alpha);
            printf("Novo h: %.8e\n\n", h);
#endif
            // usleep(300*1000);
        }
        t = t + h;
        X = vector_copy(xi, X);

        alpha = pow((eps)/(c_security*maxtal), 1.0/4.0);
        h = alpha * h;
        h = constrain(h, hmin, hmax);
        h = min(h, tf - t);

#if DEBUG > 0
        printf("Alfa:   %.8e\n", alpha);
        printf("ti:     %.8e  h: %.8e\n", t, h);
        printf("Novo h: %.8e\n\n", h);
#endif

        fprintf(out, "%.15e", t);
        for (int i = 0; i < xi->size; i++)
            fprintf(out, "  %.15e", vector_get(xi, i));
        fprintf(out, "\n");
    }
    fclose(out);
}
