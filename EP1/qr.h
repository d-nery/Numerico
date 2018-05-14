/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _QR_H_
#define _QR_H_

#include "matrix.h"
#include "vector.h"

vector_t* system_solve(matrix_t* R, vector_t* x, vector_t* bt);
void householder(matrix_t* A, vector_t* b);

#endif
