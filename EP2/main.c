/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 *
 * Exercicio Programa 2
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 *         Mateus Almeida Barbosa        - 9349072
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rkf.h"
#include "vector.h"
#include "error.h"
#include "matrix.h"

#define M_PI 3.141592653589793

vector_t* F1(double t, vector_t* X, vector_t* res) {
    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    for (int i = 0; i < X->size; i++)
        vector_set(res, i, 1.0 + pow(vector_get(X, i) - t, 2.0));

    return res;
}

vector_t* F2(double t, vector_t* X, vector_t* res) {
    (void)t; // Avoids unused warning

    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    matrix_t* A = matrix_create(4, 4);
    matrix_set(A, 0, 0, -2);
    matrix_set(A, 0, 1, -1);
    matrix_set(A, 0, 2, -1);
    matrix_set(A, 0, 3, -2);
    matrix_set(A, 1, 0,  1);
    matrix_set(A, 1, 1, -2);
    matrix_set(A, 1, 2,  2);
    matrix_set(A, 1, 3, -1);
    matrix_set(A, 2, 0, -1);
    matrix_set(A, 2, 1, -2);
    matrix_set(A, 2, 2, -2);
    matrix_set(A, 2, 3, -1);
    matrix_set(A, 3, 0,  2);
    matrix_set(A, 3, 1, -1);
    matrix_set(A, 3, 2,  1);
    matrix_set(A, 3, 3, -2);

    res = vector_mult_matrix(A, X, res);

    matrix_free(A);

    return res;
}

vector_t* F3(double t, vector_t* X, vector_t* res) {
    (void)t; // Avoids unused warning

    if (res != VEC_NULL)
        vector_free(res);
    res = vector_create(X->size);

    matrix_t* A = matrix_create(X->size, X->size);
    for (int i = 0; i < A->l; i++) {
        for (int j = 0; j < A->c; j++) {
            if (i == j)
                matrix_set(A, i, j, -2.0);
            else if ((i < A->l - 1 && j == i + 1) || (j < A->c - 1 && i == j + 1))
                matrix_set(A, i, j, 1);
        }
    }

    res = vector_mult_matrix(A, X, res);

    matrix_free(A);

    return res;
}

void x1t(double t) {
    return (double)(t + 1)/(double)(1 - t);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s caso\n", argv[0]);
        printf("    1: Unidimensional:            Teste 1\n");
        printf("    2: Multidimensional:          Teste 2\n");
        printf("    3: Multidimensional Variavel: Teste 3\n");
        printf("        Uso: %s 3 [m] (Padrao: 7)\n", argv[0]);
        printf("    4: Circuito de Chua\n");
        printf("        Uso: %s 4 [arquivo_dados] (Padrao: chua.txt)\n\n", argv[0]);
        printf("Com o makefile use 'make' para executar todos os casos e gerar os graficos\n");
        exit(0);
    }

    // Dados para o calculo em RKF
    vector_t* X0 = NULL;
    double t0    = 0.0;
    double tf    = 0.0;
    double eps   = 0.0;
    double h     = 0.0;

    int m_caso3  = 7;

	clock_t beg = clock();

    switch(argv[1][0]) {
        case '1':
            X0 = vector_create(1);
            vector_set(X0, 0,  -18.95);
            t0  = 1.05;
            tf  = 3.0;
            eps = 1e-5;
            h   = 0.1;

            rkf45_solve(X0, t0, tf, eps, h, F1, "out1.txt");
            rkf45_error(x1t, "out1.txt");
            break;

        case '2':
            X0 = vector_create(4);
            vector_set(X0, 0,  1);
            vector_set(X0, 1,  1);
            vector_set(X0, 2,  1);
            vector_set(X0, 3, -1);

            t0  = 0.0;
            tf  = 2.0;
            eps = 1e-5;
            h   = 0.1;

            rkf45_solve(X0, t0, tf, eps, h, F2, "out2.txt");
            break;

        case '3':
            if (argc > 2)
                if ((m_caso3 = atoi(argv[2])) <= 0)
                    error(ERR_NEGATIVE, "main");

            X0 = vector_create(m_caso3);
            for (int i = 0; i < X0->size; i++)
                vector_set(X0, i, sin(M_PI*((double)(i+1)/((double)m_caso3 + 1.0))) + sin(m_caso3*M_PI*((double)(i + 1)/((double)m_caso3 + 1.0))));

            t0  = 0.0;
            tf  = 2.0;
            eps = 1e-5;
            h   = 0.1;

            rkf45_solve(X0, t0, tf, eps, h, F3, "out3.txt");
            break;

        case '4':
            printf("Chua ainda nao implementado\n");
            break;

        default:
            printf("Erro! Escolha um caso valido\n");
    }
	printf("Finished! Time: %.8lfs\n\n", (double)(clock() - beg)/CLOCKS_PER_SEC);
    return 0;
}
