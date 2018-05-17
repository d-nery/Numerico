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

/**
 * newton(vector_t* ()(vector_t*), matrix_t* ()(vector_t*), vector_t*)
 * Acha raiz das funções de F(x) pelo método de Newton
 * Utiliza decomposição LU para resolver sistemas Jx = -F multiplas vezes
 *
 * Parametros:
 *  F    -> Função F(x)
 *  J    -> Função da matriz jacabiana J(x)
 *  x    -> Vetor com chutes iniciais
 *
 * Retorna
 *  x -> Modifica vetor dos valores iniciais com o resultado e o retorna
 */
vector_t* newton(vector_t* F(vector_t*), matrix_t* J(vector_t*), vector_t* x);

#endif
