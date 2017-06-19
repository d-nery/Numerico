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

const double h_coeff[]  = { 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0 };
const double k2_coeff[] = { 1.0/4.0 };
const double k3_coeff[] = { 3.0/32.0, 9.0/32.0 };
const double k4_coeff[] = { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0 };
const double k5_coeff[] = { 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0 };
const double k6_coeff[] = { -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0 };

const double xi_coeff[]  = { 25.0/216.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0 };
const double xib_coeff[] = { 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 };

#endif
