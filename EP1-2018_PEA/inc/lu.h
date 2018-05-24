/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _LU_H_
#define _LU_H_

#include "matrix.h"
#include "vector.h"

/**
 * lu()
 * Realiza a decomposição LU da matriz A, modificando-a
 *
 * Retorna o vetor p de permutações
 */
vector_t* lu(matrix_t* A, vector_t* p);


/**
 * lu_solve()
 * Resolve o sistem Ax = b
 * A e p devems er obtidos previamente a partir
 * da decomposição LU
 *
 * Retorna o vetor x
 */
vector_t* lu_solve(matrix_t* A, vector_t* x, vector_t* b, vector_t* p);


// Versao paralela para grandes dados
void lu_parallel(double** A, int n);
void lu_solve_parallel(double** A, double* x, double* b, int n);

#endif
