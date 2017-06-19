/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _RKF_H_
#define _RKF_H_

#include "vector.h"
#include "matrix.h"

#define C1 10e-9
#define C2 100e-9
#define L  18e-3

#define E    1.17391304
#define Emax 8.1818

#define Ga -50.0e-3/66.0
#define Gb -9e-3/22
#define Gc 4.591e-3

#define c_security 2.0
#define hmin       0.001
#define hmax       1.0

#define g(v) (v) <= Emax              ? Gc*(v) + Emax*(Gc - Gb) + E*(Gb - Ga) : \
             (v) > -Emax && (v) <= -E ? Gb*(v) + (Gb - Ga)*E                  : \
             (v) > -E && (v) < E      ? Ga*(v)                                : \
             (v) >= E && (v) < Emax   ? Gb*(v) + (Ga - Gb)*E                  : \
                                        Gc*(v) + Emax*(Gb - Gc) + E*(Ga - Gb)

#define min(a, b)          ((a) < (b) ? (a) : (b))
#define max(a, b)          ((a) > (b) ? (a) : (b))
#define constrain(v, a, b) ((v) < (a) ? (a) : (v) > (b) ? (b) : (v))

typedef void (*xt_function)(double);

void rkf45_solve(vector_t* X0, double t0, double tf, double eps, double h, vector_t* f(double, vector_t*, vector_t*), char* name);
void rkf45_error(void xt(double), char* file);

#endif
