COMPILE  = gcc -Wall -Wextra -std=c99 -fshort-enums
OPTIMIZE = -O3 # -fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

all: map

%.o: %.c
	$(COMPILE) -c $< -o $@ $(OPTIMIZE)
	@echo

map: $(OBJECTS)
	$(COMPILE) -o map $(OBJECTS) $(LD_FLAGS)
	@echo

out1.txt: map
	@echo "Executando caso 1"
	@./map 1

out2.txt: map
	@echo "Executando caso 2"
	@./map 2

out3.txt: map
	@echo "Executando caso 3, com m = 7"
	@./map 3 7

run: map
	@echo "Executando caso 1"
	@./map 1
	@echo "Executando caso 2"
	@./map 2
	@echo "Executando caso 3, com m = 7"
	@./map 3 7

clean:
	rm -rf $(OBJECTS) map *.txt

PHONY: all clean plot run
