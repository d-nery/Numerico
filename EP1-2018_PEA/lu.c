#include <stdio.h>
#include <math.h>

#include "log.h"
#include "lu.h"
#include "matrix.h"
#include "vector.h"
#include "error.h"

void lu(matrix_t* A) {
    if (A == MAT_NULL)
        error(ERR_NULL, "lu");

    if (A->l != A->c)
        error(ERR_SIZE, "lu");

    // if (L != MAT_NULL)
    //     matrix_free(L);
    // L = matrix_create(A->l, A->c);

    // if (U != MAT_NULL)
    //     matrix_free(U);
    // U = matrix_create(A->l, A->c);

    const int n = A->l;
    double sum = 0;

    vector_t* p = vector_create(n);

    for (int k = 0; k < n; k++) {
        log_info("K = %d", k);
        // Passo 1
        for (int i = k; i < n; i++) {
            log_info("I = %d", i);

            sum = 0;
            for (int j = 0; j < k - 1; j++) {
                sum += matrix_get(A, i, j) * matrix_get(A, j, k);
            }
            log_trace("sum -> %f", sum);
            matrix_set(A, i, k, matrix_get(A, i, k) - sum);

            print_matrix(A);
        }

        // Determinar l
        double _max = 0, _old_max = 0;
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
            for (int i = 0; j < k - 1; i++) {
                sum += matrix_get(A, k, i) * matrix_get(A, i, j);
            }
            matrix_set(A, k, j, matrix_get(A, k, j) - sum);
            matrix_set(A, j, k, matrix_get(A, j, k)/matrix_get(A, k, k));
        }
    }
}
