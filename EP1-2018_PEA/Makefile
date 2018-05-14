COMPILE  = gcc -Wall -Wextra -std=c11 -fshort-enums
OPTIMIZE = -O3 # -fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

CFLAGS = -DDEBUG=5

all: main

%.o: %.c
	$(COMPILE) -c $< -o $@ $(OPTIMIZE) $(CFLAGS)
	@echo

main: $(OBJECTS)
	$(COMPILE) -o main $(OBJECTS) $(LD_FLAGS)
	@echo

clean:
	rm -rf $(OBJECTS) main main.exe

.PHONY: all clean
