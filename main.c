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

#define len(array) ((&array)[1] - array)

int main(int argc, char* argv[]) {
	matrix_t* A = matrix_create_from_file("data/6259_Completa_D_Matriz.txt");
	printf("\n");

	matrix_t* B = matrix_plus(A, A);

	return 0;
}
