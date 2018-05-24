/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _TESTES_H_
#define _TESTES_H_

#include "matrix.h"
#include "vector.h"

/**
 * Testes iniciais
 * Funções passadas para o metodo de newton
 */

/**
 * F1
 * F e J do teste inicial 1
 */
vector_t* F1(vector_t* x);
matrix_t* JF1(vector_t* x);

/**
 * F2
 * F e J do teste inicial 2
 */
vector_t* F2(vector_t* x);
matrix_t* JF2(vector_t* x);

/**
 * F3
 * F e J do teste inicial 3
 */
vector_t* F3(vector_t* x);
matrix_t* JF3(vector_t* x);

#endif
