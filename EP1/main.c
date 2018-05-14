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
#include "qr.h"

#define min(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: %s [matrix_size_file]\n", argv[0]);
		exit(0);
	}

	char file_m[64];
	char file_v[64];
	sprintf(file_m, "data/%s_Completa_D_Matriz.txt", argv[1]);
	sprintf(file_v, "data/%s_Completa_D_VetorB.txt", argv[1]);

	matrix_t* A = matrix_create_from_file(file_m);
	vector_t* b = vector_create_from_file(file_v, A->l);

	printf("Starting computing...\n");
	clock_t beg = clock();
	householder(A, b);

	vector_t* x = NULL;
	x = system_solve(A, x, b);

	printf("Finished! Time: %.8lfs\n", (double)(clock() - beg)/CLOCKS_PER_SEC);
	sprintf(file_v, "out/%s_X.txt", argv[1]);
	output_vector(x, file_v);

	return 0;
}
