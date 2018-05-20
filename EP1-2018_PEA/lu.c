/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "log.h"
#include "lu.h"
#include "matrix.h"
#include "vector.h"
#include "error.h"

vector_t* lu(matrix_t* A, vector_t* p) {
    log_trace("Iniciando decomposição LU");

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
        // Passo 1
        for (int i = k; i < n; i++) {
            sum = 0;
            for (int j = 0; j < k; j++) {
                sum += matrix_get(A, i, j) * matrix_get(A, j, k);
            }
            matrix_set(A, i, k, matrix_get(A, i, k) - sum);
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

        // Trocar linhas
        if (k != (int)vector_get(p, k)) {
            matrix_swap_lines(A, k, (int)vector_get(p, k));
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

        if (k % 100 == 0)
            log_trace("LU k = %d/%d (%d%%)", k, n, (int)((100.0 * k) / n));
    }

    log_trace("LU k = %d/%d (%d%%)", n, n, 100);
    log_trace("Concluido");

    return p;
}

vector_t* lu_solve(matrix_t* A, vector_t* x, vector_t* b, vector_t* p) {
    log_trace("Iniciando resolução do sistema");

    if (A == MAT_NULL || b == VEC_NULL || p == VEC_NULL)
        error(ERR_NULL, "lu_solve");

    if (x == VEC_NULL) {
        x = vector_create(b->size);
    }

    double temp;

    // Permuta b
    for (int i = 0; i < p->size; i++) {
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

    log_trace("Concluido");

    return x;
}

void lu_parallel(double** A, int n) {
    // Essa versão não faz pivotação, não é problema
    // para as redes, mas usar com cuidado
    log_trace("Iniciando decomposição LU");

    int i, k;

    for (k = 0; k < n - 1; k++) {
        #pragma omp parallel for
        for (i = k + 1; i < n; i++) {
            A[i][k] /= A[k][k];
            if (isnan(A[i][k]) || isinf(A[i][k]))
                log_warn("NaN or Inf found [%d, %d]", i, k);
        }

        #pragma omp parallel for shared(A, n, k) private(i) schedule(static, 100) num_threads(8)
        for (i = k + 1; i < n; i++) {
            int j;
            for (j = k + 1; j < n; j++) {
                A[j][i] -= A[k][i] * A[j][k];
            }
        }

        if (k % 100 == 0)
            log_trace("LU k = %d/%d (%d%%)", k, n, (int)((100.0 * k) / n));
    }

    log_trace("LU k = %d/%d (%d%%)", k, n, (int)((100.0 * k+1) / n));

    log_trace("Concluido");
}


void lu_solve_parallel(double** A, double* x, double* b, int n) {
    log_trace("Iniciando resolução do sistema");

    for (int i = 0; i < n; i++) {
        x[i] = b[i];

        for (int j = 0; j < i; j++) {
            x[i] -= A[i][j] * x[j];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }

        x[i] /= A[i][i];
    }

    log_trace("Concluido");
}
