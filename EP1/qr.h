/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */


#include "matrix.h"
#include "vector.h"

void r_solve(matrix_t* R, vector_t* x, vector_t* bt);
void qr_factorize(matrix_t* A, vector_t* b);
