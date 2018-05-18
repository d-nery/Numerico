/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _REDE_H_
#define _REDE_H_

#include "matrix.h"
#include "vector.h"

/**
 * Funções para calculo das Redes propostas
 * prepara_rede -> Le os arquivos relacionados a rede n
 *                 e inicializa o que for necessario
 * F_rede -> retorna o F(x) para a rede preparada
 * J_rede -> retorna o J(x) para a rede preparada
 * finaliza_rede -> Chamada apos o achar as raizes por newton
 */

vector_t* F_rede(vector_t* x);
matrix_t* J_rede(vector_t* x);

vector_t* prepara_rede(const int n);
void finaliza_rede(vector_t* x);

#endif
