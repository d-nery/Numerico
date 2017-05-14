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
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "matrix.h"
#include "vector.h"
#include "error.h"


int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: %s [matrix_size_file]\n", argv[0]);
		exit(0);
	}

	char file_m[64];
	char file_v[64];
	sprintf(file_m, "data/%s_Completa_D_Matriz.txt", argv[1]);
	sprintf(file_v, "data/%s_Completa_D_VetorB.txt", argv[1]);

	// matrix_t* A = matrix_create_from_file(file_m);
	// vector_t* b = vector_create_from_file(file_v, A->l);

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	matrix_t* A = matrix_create(n, m);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			A->data[i][j] = abs(i - j) <= 4 ? 1.0/(i + j -1) : 0;

	vector_t* b = vector_create(n);
	for (int i = 0; i < n; i++)
		b->data[i] = 1;

	printf("Starting computing...\n");
	clock_t beg = clock();
	for (int c = 0; c < A->c-1; c++) {
		vector_t* B = vector_create(A->l - c);
		for (int j = 0; j < B->size; j++)
			vector_set(B, j, matrix_get(A, j + c, c));

		vector_t* e  = vector_create(B->size);
		vector_set(e, 0, 1.0);


		vector_t* hb = vector_create(B->size);
		vector_t* w  = NULL;
		vector_t* temp_mult = NULL;
		temp_mult = vector_mult_scalar(vector_norm(B), e, temp_mult);
		w = vector_add(B, temp_mult, w);

		for (int k = 0; k < B->size; k++) {
			for (int j = 0; j < B->size; j++)
				vector_set(B, j, matrix_get(A, j + c, k + c));

			// H*b = B - 2 * (w * B)/(w * w) * w
			temp_mult = vector_mult_scalar(2 * vector_multiply(w, B)/vector_multiply(w, w), w, temp_mult);
			hb = vector_subtract(B, temp_mult, hb);

			for (int i = c; i < A->c; i++)
				matrix_set(A, i, k + c, vector_get(hb, i - c));

		}
		// Vetor b
		for (int j = 0; j < B->size; j++)
			vector_set(B, j, vector_get(b, j + c));

		temp_mult = vector_mult_scalar(2 * vector_multiply(w, B)/vector_multiply(w, w), w, temp_mult);
		hb = vector_subtract(B, temp_mult, hb);
		for (int i = c; i < A->c; i++)
			vector_set(b, i, vector_get(hb, i - c));

		vector_free(B);
		vector_free(w);
		vector_free(e);
		vector_free(hb);
		vector_free(temp_mult);
	}

	vector_t* x = vector_create(b->size);

	for (int i = A->l-1; i >= 0; i--) {
		double val = 0;
		for (int j = i+1; j < A->c; j++) {
			val += vector_get(x, j) * matrix_get(A, i, j);
		}
		vector_set(x, i, (vector_get(b, i) - val)/matrix_get(A, i, i));
	}
	printf("Finished! Time: %.8lfs\n", (double)(clock() - beg)/CLOCKS_PER_SEC);
	sprintf(file_v, "out/%s_X.txt", argv[1]);
	output_vector(x, file_v);

	return 0;
}
