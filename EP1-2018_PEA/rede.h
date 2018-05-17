#ifndef _REDE_H_
#define _REDE_H_

#include "matrix.h"
#include "vector.h"

vector_t* F_rede(vector_t* x);
matrix_t* J_rede(vector_t* x);

vector_t* prepara_rede(const int n);

#endif
