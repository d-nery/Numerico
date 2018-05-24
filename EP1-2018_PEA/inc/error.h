/**
 * Escola Politécnica da USP
 * MAP3121 - Metodos Numericos e Aplicacoes
 * PEA3301 - Introdução aos Sistemas de Potencia
 *
 * Exercicio Programa 1
 *
 * Alunos: Daniel Nery Silva de Oliveira - 9349051
 */

#ifndef _ERROR_H_
#define _ERROR_H_

typedef enum {
    ERR_UNKNOWN = 0,
    ERR_NEGATIVE,
    ERR_MEMORY,
    ERR_OPEN_FILE,
    ERR_READ_FILE,
    ERR_SIZE,
    ERR_OOB,
    ERR_NULL
} err_t;

#define error(err, func) __err(err, func, __FILE__, __LINE__)
void __err(err_t err, char* func, char* file, int line);

#endif
