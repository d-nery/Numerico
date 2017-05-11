#include <stdio.h>
#include <stdlib.h>

#include "error.h"

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

struct {
	char** list;
	int len;
} errors = { err_list, 7 };

void __err(err_t err, char* func, char* file, int line) {
	if (err >= errors.len || err < 0)
		err = 0;

	fprintf(stderr, "[ERRO] %s:%d: Na funcao %s(): %s\n", file, line, func, errors.list[err]);
	exit(0);
}
