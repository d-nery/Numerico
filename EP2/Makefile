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

plot: out1.txt out2.txt out3.txt
	@echo "Gerando graficos..."
	@gnuplot -e "set term png; set output 'out1.png'; plot 'out1.txt', x+1/(1-x)"
	@gnuplot -e "set term png; set output 'out2.png'; plot 'out2.txt' using 1:2, 'out2.txt' using 1:3, 'out2.txt' using 1:4, 'out2.txt' using 1:5, exp(-x)*sin(x)+exp(-3*x)*cos(3*x), exp(-x)*cos(x)+exp(-3*x)*sin(3*x), -exp(-x)*sin(x)+exp(-3*x)*cos(3*x), -exp(-x)*cos(x)+exp(-3*x)*sin(3*x)"
	@gnuplot -e "set term png; set output 'out3.png'; unset key; \
		plot 'out3.txt' using 1:2, \
		'out3.txt' using 1:3, \
		'out3.txt' using 1:4, \
		'out3.txt' using 1:5, \
		'out3.txt' using 1:6, \
		'out3.txt' using 1:7, \
		'out3.txt' using 1:8, \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*0/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*0/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*1/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*1/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*2/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*2/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*3/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*3/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*4/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*4/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*5/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*5/8), \
		exp(-(2*(1-cos(pi/8)))*x)*sin(pi*6/8)+exp(-(2*(1-cos(7*pi/8)))*x)*sin(7*pi*6/8)"

clean:
	rm -rf $(OBJECTS) map *.txt

PHONY: all clean plot run
