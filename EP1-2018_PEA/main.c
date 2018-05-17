#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "log.h"
#include "matrix.h"
#include "vector.h"
#include "error.h"
#include "lu.h"
#include "newton.h"
#include "testes.h"
#include "rede.h"

#ifndef DEBUG
#define DEBUG 0
#endif

void usage(char* argv[]);

int main(int argc, char* argv[]) {
    log_set_level(LOG_FATAL - DEBUG);

    if (argc < 2) {
        usage(argv);
        exit(0);
    }

    if ((argv[1][0] == '4' || argv[1][0] == '3') && argc < 3) {
        usage(argv);
        exit(0);
    }

    vector_t* x = VEC_NULL;

    switch (argv[1][0]) {
        case '0':
        {
            log_info("Teste F0:");

            matrix_t* A = matrix_create(3, 3);

            matrix_set(A, 0, 0, 1);
            matrix_set(A, 0, 1, -2);
            matrix_set(A, 0, 2, 3);

            matrix_set(A, 1, 0, 2);
            matrix_set(A, 1, 1, 1);
            matrix_set(A, 1, 2, 1);

            matrix_set(A, 2, 0, -3);
            matrix_set(A, 2, 1, 2);
            matrix_set(A, 2, 2, -2);

            vector_t* b = vector_create(3);

            vector_set(b, 0, 7);
            vector_set(b, 1, 4);
            vector_set(b, 2, -10);

            vector_t* p = VEC_NULL;

            p = lu(A, p);
            print_matrix(A);

            x = lu_solve(A, x, b, p);
            print_vector(x);

            break;
        }

        case '1':
        {
            log_info("Teste F1:");

            x = vector_create(2);

            vector_set(x, 0, 0);
            vector_set(x, 1, 0);

            x = newton(F1, JF1, x);

            print_vector(x);

            vector_free(x);
            break;
        }

        case '2':
        {
            log_info("Teste F2:");

            x = vector_create(4);

            vector_set(x, 0, 1);
            vector_set(x, 1, 1);
            vector_set(x, 2, 1);
            vector_set(x, 3, 1);

            x = newton(F2, JF2, x);

            print_vector(x);

            vector_free(x);
            break;
        }

        case '3':
        {
            log_info("Teste F3:");

            int n = atoi(argv[2]);

            x = vector_create(n - 1);

            x = newton(F3, JF3, x);

            print_vector(x);

            vector_free(x);
            break;
        }

        case '4':
        {
            log_info("Teste F4:");

            int n = atoi(argv[2]);

            x = prepara_rede(n);
            newton(F_rede, J_rede, x);
            // print_vector(x);as
            log_info("Terminado!");

            vector_free(x);

            break;
        }

        default:
            usage(argv);
    }

    return 0;
}

void usage(char* argv[]) {
    printf("Usage:\n");
    printf("    %s numero_teste [n | caso]\n\n", argv[0]);
    printf("    numero_teste:\n");
    printf("     1: Teste 1\n");
    printf("     2: Teste 2\n");
    printf("     3 n: Teste 3\n");
    printf("     4 caso: Redes\n\n");
    printf("    n: sistema n-1 X n-1\n\n");
    printf("    caso:\n");
    printf("     1: Stevenson\n");
    printf("     2: Reticulada\n");
    printf("     3: Distribuicao Primaria\n");
    printf("     4: Distribuicao Primaria Secundária\n");
}
