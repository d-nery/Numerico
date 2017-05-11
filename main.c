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
#include "error.h"

int main(int argc, char* argv[]) {
	matrix_t* A = matrix_create_from_file("data/4_Completa_D_Matriz.txt");
	printf("\n");

	print_matrix(matrix_create_identity(4));

	return 0;
}
