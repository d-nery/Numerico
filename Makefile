<<<<<<< HEAD
COMPILE  = gcc -Wall -Wextra -std=c99
OPTIMIZE = -Ofast -fopenmp
=======
COMPILE  = gcc -Wall -Wextra -std=c99 -fast
>>>>>>> 24bf126834d9e7bbd6e1638bbda2a65dd4b2ea5b
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

all: map

%.o: %.c
	$(COMPILE) $(OPTIMIZE) -c $< -o $@
	@echo

map: $(OBJECTS)
	$(COMPILE) -o map $(LD_FLAGS) $(OBJECTS)
	@echo

clean:
	rm -rf $(OBJECTS) map

PHONY: all clean
