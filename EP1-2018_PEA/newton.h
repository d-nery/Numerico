/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _NEWTON_H_
#define _NEWTON_H_

#include "vector.h"
#include "matrix.h"

typedef vector_t* (*func)(vector_t*);
typedef matrix_t* (*jacobian_func)(vector_t*);

vector_t* newton(func F, jacobian_func J, vector_t* x);

#endif
