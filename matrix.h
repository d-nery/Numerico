/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct {
	float** data;
	int l, c;
} matrix_t;

matrix_t* matrix_create(int m, int n);

#endif
