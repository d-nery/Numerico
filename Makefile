COMPILE  = gcc -Wall -Wextra -std=c99 -fshort-enums
OPTIMIZE = -Ofast -fopenmp
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
