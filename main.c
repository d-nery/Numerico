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
#include "error.h"


int main(int argc, char* argv[]) {
	matrix_t* A = matrix_create_from_file("data/4_Completa_D_Matriz.txt");
	vector_t* B = vector_create_from_file("data/4_Completa_D_VetorB.txt", A->l);

	print_matrix(A);
	printf("\n");
	print_vector(B);
	printf("\n");

	return 0;
}
