# Escola Politécnica da USP
# MAP3121 - Metodos Numericos e Aplicacoes
# PEA3301 - Introdução aos Sistemas de Potencia
#
# Exercicio Programa 1
#
# Alunos: Daniel Nery Silva de Oliveira - 9349051

COMPILE  = gcc -Wall -Wextra -std=c11 -fshort-enums
OPTIMIZE = -O2 -fopenacc # -fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

# DEBUG: 0 -> FATAL, 1 -> ERROR, 2-> WARNING, 3 -> INFO, 4-> DEBUG, 5 -> TRACE
CFLAGS = -DDEBUG=5
ifdef SYSTEMROOT # Windows, no color prompt
CFLAGS +=
else
CFLAGS += -DLOG_USE_COLOR
endif

all: main

%.o: %.c
	$(COMPILE) -c $< -o $@ $(OPTIMIZE) $(CFLAGS) -I"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.2\include"
	@echo

main: $(OBJECTS)
	$(COMPILE) -o main $(OBJECTS) $(LD_FLAGS) $(OPTIMIZE) -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.2\lib\x64" "-Wl,-rpath,C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.2\lib\x64"
	@echo

clean:
	rm -rf $(OBJECTS) main main.exe

.PHONY: all clean
