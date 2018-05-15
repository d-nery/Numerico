#include <math.h>

#include "newton.h"
#include "log.h"
#include "error.h"
#include "vector.h"
#include "matrix.h"
#include "lu.h"

#define EPS 1e-15

vector_t* newton(func F, jacobian_func J, vector_t* x) {
    if (x == VEC_NULL)
        error(ERR_NULL, "newton_solve");

    vector_t* p = VEC_NULL;
    vector_t* c = VEC_NULL;

    vector_t* Fx = VEC_NULL;
    matrix_t* Jx = MAT_NULL;

    int it = 0;

    for (;;) {
        Fx = F(x);
        Jx = J(x);

        Fx = vector_mult_scalar(-1, Fx, Fx);

        p = lu(Jx, p);
        c = lu_solve(Jx, c, Fx, p);

        x = vector_add(x, c, x);

        if (vector_norm(c) < EPS)
            break;

        it++;
    }

    log_trace("Newton terminado: %d iteracoes", it);

    return x;
}
