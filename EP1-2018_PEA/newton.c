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
        log_trace("iteração %d", it+1);
        log_trace("Calculando F");
        Fx = F(x);
        log_trace("Concluido");
        log_trace("Calculando J");
        Jx = J(x);
        log_trace("Concluido");

        Fx = vector_mult_scalar(-1, Fx, Fx);

        log_trace("Resolvendo Sistema");
        p = lu(Jx, p);
        c = lu_solve(Jx, c, Fx, p);
        log_trace("Concluido");

        x = vector_add(x, c, x);

        it++;

        vector_free(Fx);
        matrix_free(Jx);

        log_trace("Checando convergencia");
        if (vector_norm(c) < EPS)
            break;
    }

    log_debug("Newton terminado: %d iteracoes", it);

    return x;
}
