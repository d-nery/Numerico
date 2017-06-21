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

#define c_security 2.0
#define hmin       1e-9
#define hmax       1.0

#define min(a, b)          ((a) < (b) ? (a) : (b))
#define max(a, b)          ((a) > (b) ? (a) : (b))
#define constrain(v, a, b) ((v) < (a) ? (a) : (v) > (b) ? (b) : (v))

/**
 * rkf45_solve(vector_t*, double, double, double, double, vector_t* ()(double, vector_t*, vector_t*), char*)
 * Soluciona equacao diferencial pelo metodo RKF45
 *
 * Parametros:
 *  X0   -> vetor com os valores iniciais
 *  t0   -> tempo inicial
 *  tf   -> tempo final
 *  eps  -> precisão das iterações
 *  h    -> h inicial
 *  f    -> funcao F(t, X)
 *  xt   -> funcao X(t), com resultado real da equacao, para calculo do erro
 *  name -> nome do arquivo de saida dos dados
 */
void rkf45_solve(vector_t* X0, double t0, double tf, double eps, double h, vector_t* f(double, vector_t*, vector_t*), char* name);


#endif
