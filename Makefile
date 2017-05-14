COMPILE  = gcc -Wall -Wextra -std=c99 -fshort-enums -g
OPTIMIZE = -Og -fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

all: map

%.o: %.c
	$(COMPILE) -c $< -o $@ $(OPTIMIZE)
	@echo

map: $(OBJECTS)
	$(COMPILE) -o map $(OBJECTS) $(LD_FLAGS)
	@echo

clean:
	rm -rf $(OBJECTS) map

PHONY: all clean
