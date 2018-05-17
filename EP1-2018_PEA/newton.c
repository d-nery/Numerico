/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#include <math.h>

#include "newton.h"
#include "log.h"
#include "error.h"
#include "vector.h"
#include "matrix.h"
#include "lu.h"

#define EPS 1e-3

vector_t* newton(vector_t* F(vector_t*), matrix_t* J(vector_t*), vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "newton_solve");

    vector_t* p = VEC_NULL;
    vector_t* c = VEC_NULL;

    vector_t* Fx = VEC_NULL;
    matrix_t* Jx = MAT_NULL;

    int it = 0;

    for (;;) {
        // log_info("IT: %d", it+1);
        Jx = J(x);
        Fx = F(x);

        Fx = vector_mult_scalar(-1, Fx, Fx);

        p = lu(Jx, p);
        c = lu_solve(Jx, c, Fx, p);
        // log_info("C");
        // print_vector(c);

        x = vector_add(x, c, x);

        if (vector_norm(c) < EPS)
            break;

        // print_vector(x);

        vector_free(Fx);
        matrix_free(Jx);

        it++;
    }

    log_trace("Newton terminado: %d iteracoes", it);

    return x;
}
