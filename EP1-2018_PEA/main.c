#include <stdio.h>
#include <stdlib.h>

#include "log.h"
#include "matrix.h"
#include "lu.h"

#ifndef DEBUG
#define DEBUG 0
#endif

int main(int argc, char* argv[]) {
    log_set_level(LOG_FATAL - DEBUG);

    matrix_t* A = matrix_create(3,3);
    matrix_t* L = matrix_create(3,3);
    matrix_t* U = matrix_create(3,3);

    matrix_t* R;

    matrix_set(A, 0, 0, 1);
    matrix_set(A, 0, 1, 2);
    matrix_set(A, 0, 2, 3);
    matrix_set(A, 1, 0, 4);
    matrix_set(A, 1, 1, 5);
    matrix_set(A, 1, 2, 6);
    matrix_set(A, 2, 0, 7);
    matrix_set(A, 2, 1, 8);
    matrix_set(A, 2, 2, 9);

    print_matrix(A);
    lu(A);

    L->data[0][0] = 1;
    L->data[1][0] = A->data[1][0];
    L->data[1][1] = 1;
    L->data[2][0] = A->data[2][0];
    L->data[2][1] = A->data[2][1];
    L->data[2][2] = 1;

    U->data[0][0] = A->data[0][0];
    U->data[0][1] = A->data[0][1];
    U->data[0][2] = A->data[0][2];
    U->data[1][1] = A->data[1][1];
    U->data[1][2] = A->data[1][2];
    U->data[2][2] = A->data[2][2];

    log_info("Imprimindo resultado: ");
    print_matrix(L);
    print_matrix(U);
    R = matrix_multiply(L, U, R);
    print_matrix(R);

    return 0;
}
