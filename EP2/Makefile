COMPILE  = gcc -Wall -Wextra -std=c99 -fshort-enums
OPTIMIZE = -O3 #-fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.c,%.o,$(wildcard *.c))

all: map

%.o: %.c
	$(COMPILE) -c $< -o $@ $(OPTIMIZE)
	@echo

map: $(OBJECTS)
	$(COMPILE) -o map $(OBJECTS) $(LD_FLAGS)
	@echo

out.txt: map

plot: out.txt out2.txt
	@gnuplot -e "set term png; set output 'out1.png'; plot 'out.txt', x+1/(1-x); set term x11"
	@gnuplot -e "set term png; set output 'out2.png'; plot 'out2.txt' using 1:2, 'out2.txt' using 1:3, 'out2.txt' using 1:4, 'out2.txt' using 1:5, exp(-x)*sin(x)+exp(-3*x)*cos(3*x), exp(-x)*cos(x)+exp(-3*x)*sin(3*x), -exp(-x)*sin(x)+exp(-3*x)*cos(3*x), -exp(-x)*cos(x)+exp(-3*x)*sin(3*x); set term x11"

clean:
	rm -rf $(OBJECTS) map.exe

PHONY: all clean plot
