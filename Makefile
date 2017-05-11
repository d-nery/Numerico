COMPILE  = gcc -Wall -Wextra -std=c99 -fast
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

all: map

%.o: %.c
	$(COMPILE) -c $< -o $@
	@echo

map: $(OBJECTS)
	$(COMPILE) -o map $(LD_FLAGS) $(OBJECTS)
	@echo

clean:
	rm -rf $(OBJECTS) map.exe

PHONY: all clean
