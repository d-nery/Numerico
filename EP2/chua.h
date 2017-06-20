/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _CHUA_H_
#define _CHUA_H_

#define C1 10e-9
#define C2 100e-9
#define L  18e-3
#define R  1000.0

#define E    1.17391304
#define Emax 8.1818

#define Ga -50.0e-3/66.0
#define Gb -9e-3/22
#define Gc 4.591e-3

// #define g(v) (v) <= -Emax              ? (double)(Gc*(v) + Emax*(Gc - Gb) + E*(Gb - Ga)) : \
//              (v) >  -Emax && (v) <= -E ? (double)(Gb*(v) + (Gb - Ga)*E)                  : \
//              (v) >  -E && (v) < E      ? (double)(Ga*(v))                                : \
//              (v) >=  E && (v) < Emax   ? (double)(Gb*(v) + (Ga - Gb)*E)                  : \
//                                          (double)(Gc*(v) + Emax*(Gb - Gc) + E*(Ga - Gb))

inline double g(double v) {
    return
        (v) <= -Emax              ? (Gc*(v) + Emax*(Gc - Gb) + E*(Gb - Ga)) :
        (v) >  -Emax && (v) <= -E ? (Gb*(v) + (Gb - Ga)*E)                  :
        (v) >  -E && (v) < E      ? (Ga*(v))                                :
        (v) >=  E && (v) < Emax   ? (Gb*(v) + (Ga - Gb)*E)                  :
                                    (Gc*(v) + Emax*(Gb - Gc) + E*(Ga - Gb));
}

#endif
