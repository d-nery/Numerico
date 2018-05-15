#ifndef _LU_H_
#define _LU_H_

#include "matrix.h"
#include "vector.h"

vector_t* lu(matrix_t* A, vector_t* p);
vector_t* lu_solve(matrix_t* A, vector_t* x, vector_t* b, vector_t* p);

#endif
