COMPILE  = g++ -Wall -Wextra -std=c++1z
OPTIMIZE = -Ofast -fopenmp
LD_FLAGS = -lm
OBJECTS  = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

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
