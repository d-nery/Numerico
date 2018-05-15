#ifndef _NEWTON_H_
#define _NEWTON_H_

#include "vector.h"
#include "matrix.h"

typedef vector_t* (*func)(vector_t*);
typedef matrix_t* (*jacobian_func)(vector_t*);

vector_t* newton(func F, jacobian_func J, vector_t* x);

#endif
