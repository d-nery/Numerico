#include <stdio.h>
#include <math.h>

#include "log.h"
#include "lu.h"
#include "matrix.h"
#include "vector.h"
#include "error.h"

vector_t* lu(matrix_t* A, vector_t* p) {
    if (A == MAT_NULL)
        error(ERR_NULL, "lu");

    if (A->l != A->c)
        error(ERR_SIZE, "lu");

    if (p != VEC_NULL)
        vector_free(p);
    p = vector_create(A->l);

    const int n = A->l;
    double sum = 0;

    for (int k = 0; k < n; k++) {
        log_info("K = %d", k);
        // Passo 1
        for (int i = k; i < n; i++) {
            log_info("I = %d", i);

            sum = 0;
            for (int j = 0; j < k; j++) {
                sum += matrix_get(A, i, j) * matrix_get(A, j, k);
            }
            log_trace("sum -> %f", sum);
            matrix_set(A, i, k, matrix_get(A, i, k) - sum);

            print_matrix(A);
        }

        // Determinar l
        double _max = 0, _old_max = -1;
        for (int i = k; i < n; i++) {
            _max = fmax(_max, fabs(matrix_get(A, i, k)));
            if (_max != _old_max) {
                vector_set(p, k, i);
            }
            _old_max = _max;
        }

        print_vector(p);

        // Trocar linhas
        if (k != (int)vector_get(p, k)) {
            matrix_swap_lines(A, k, (int)vector_get(p, k));
            log_info("Linhas trocadas! %d <-> %d", k, (int)vector_get(p, k));
            print_matrix(A);
        }

        // Ultimo passo
        for (int j = k + 1; j < n; j++) {
            sum = 0;
            for (int i = 0; i < k; i++) {
                sum += matrix_get(A, k, i) * matrix_get(A, i, j);
            }
            matrix_set(A, k, j, matrix_get(A, k, j) - sum);
            matrix_set(A, j, k, matrix_get(A, j, k)/matrix_get(A, k, k));
        }
    }

    return p;
}


vector_t* lu_solve(matrix_t* A, vector_t* x, vector_t* b, vector_t* p) {
    x = vector_create(b->size);
    double temp;

    // Permuta b
    for (int i = p->size-1; i >= 0; i--) {
        int ind = (int)vector_get(p, i);
        if (i != ind) {
            temp = vector_get(b, i);
            vector_set(b, i, vector_get(b, ind));
            vector_set(b, ind, temp);
        }
    }

    for (int i = 0; i < x->size; i++) {
        vector_set(x, i, vector_get(b, i));

        for (int j = 0; j < i; j++) {
            vector_set(x, i, vector_get(x, i) - matrix_get(A, i, j)*vector_get(x, j));
        }
    }

    for (int i = x->size - 1; i >= 0; i--) {
        for (int j = i + 1; j < x->size; j++) {
            vector_set(x, i, vector_get(x, i) - matrix_get(A, i, j)*vector_get(x, j));
        }

        vector_set(x, i, vector_get(x, i)/matrix_get(A, i, i));
    }

    return x;
}
