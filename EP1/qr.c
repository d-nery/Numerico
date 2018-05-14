/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <stdio.h>

#include "matrix.h"
#include "vector.h"
#include "qr.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

static int it_max;

/**
 * system_solve()
 * Resolve um sistema linear do tipo Rx = b, onde R é uma matriz
 * triangular inferior
 *
 * Retorna o vetor x
 */
vector_t* system_solve(matrix_t* R, vector_t* x, vector_t* bt) {
	// if (x != VEC_NULL)
	// 	vector_free(x);
	x = vector_create(R->c);

	for (int i = it_max; i >= 0; i--) {
		double val = 0;
		for (int j = i+1; j < R->c; j++) {
			val += vector_get(x, j) * matrix_get(R, i, j);
		}
		vector_set(x, i, (vector_get(bt, i) - val)/matrix_get(R, i, i));
	}

	return x;
}

void householder(matrix_t* A, vector_t* b) {
	it_max = min(A->l-1, A->c);
	for (int c = 0; c < it_max; c++) {
		vector_t* B = vector_create(A->l - c);
		for (int j = 0; j < B->size; j++)
			vector_set(B, j, matrix_get(A, j + c, c));

		vector_t* e  = vector_create(B->size);
		vector_set(e, 0, 1.0);

		vector_t* hb = vector_create(B->size);
		vector_t* w  = NULL;
		vector_t* temp_mult = NULL;
		temp_mult = vector_mult_scalar(sgn(vector_get(B, 0)) * vector_norm(B), e, temp_mult);
		w = vector_add(B, temp_mult, w);

		for (int k = c; k < A->c; k++) {
			// printf("--- Coluna %d\n", k);
			for (int j = 0; j < B->size; j++)
				vector_set(B, j, matrix_get(A, j + c, k));

			// H*b = B - 2 * (w * B)/(w * w) * w
			temp_mult = vector_mult_scalar(2 * vector_multiply(w, B)/vector_multiply(w, w), w, temp_mult);
			hb = vector_subtract(B, temp_mult, hb);

			for (int i = c; i < A->l; i++)
				matrix_set(A, i, k, vector_get(hb, i - c));
		}

		// Vetor b
		for (int j = 0; j < B->size; j++)
			vector_set(B, j, vector_get(b, j + c));

		temp_mult = vector_mult_scalar(2 * vector_multiply(w, B)/vector_multiply(w, w), w, temp_mult);
		hb = vector_subtract(B, temp_mult, hb);
		for (int i = c; i < A->l; i++)
			vector_set(b, i, vector_get(hb, i - c));

		vector_free(B);
		vector_free(w);
		vector_free(e);
		vector_free(hb);
		vector_free(temp_mult);
	}
}
