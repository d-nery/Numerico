#include <stdio.h>
#include <stdlib.h>

#include "error.h"

#define len(array) ((&array)[1] - array)

char* err_list[] = {
	"erro desconhecido",
	"argumento negativo",
	"nao foi possivel alocar memoria",
	"erro na abertura do arquivo",
	"erro na leitura do arquivo",
	"tamanhos nao condizentes de matriz/vetor",
	"out of bounds",
	"argumentos nulos"
};

void __err(err_t err, char* func, char* file, int line) {
	if (err >= len(err_list))
		err = 0;

	fprintf(stderr, "[ERRO] %s:%d: Na funcao %s(): %s\n", file, line, func, err_list[err]);
	exit(0);
}
