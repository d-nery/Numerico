#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "log.h"
#include "matrix.h"
#include "vector.h"
#include "lu.h"

#ifndef DEBUG
#define DEBUG 0
#endif

int main(int argc, char* argv[]) {
    log_set_level(LOG_FATAL - DEBUG);

    matrix_t* A = matrix_create(3,3);

    vector_t* b = vector_create(A->l);

    vector_t* p;
    vector_t* x;

    matrix_set(A, 0, 0, 1);
    matrix_set(A, 0, 1, -2);
    matrix_set(A, 0, 2, 3);

    matrix_set(A, 1, 0, 2);
    matrix_set(A, 1, 1, 1);
    matrix_set(A, 1, 2, 1);

    matrix_set(A, 2, 0, -3);
    matrix_set(A, 2, 1, 2);
    matrix_set(A, 2, 2, -2);

    vector_set(b, 0, 7);
    vector_set(b, 1, 4);
    vector_set(b, 2, -10);

    print_matrix(A);
    print_vector(b);

    p = lu(A, p);
    x = lu_solve(A, x, b, p);

    log_info("Imprimindo resultado: ");
    print_vector(x);

    return 0;
}
