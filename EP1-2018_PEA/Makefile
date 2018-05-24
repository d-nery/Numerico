# Escola Politécnica da USP
# MAP3121 - Metodos Numericos e Aplicacoes
# PEA3301 - Introdução aos Sistemas de Potencia
#
# Exercicio Programa 1
#
# Aluno: Daniel Nery Silva de Oliveira - 9349051

TARGET = main
# DEBUG: 0 -> FATAL, 1 -> ERROR, 2-> WARNING, 3 -> INFO, 4-> DEBUG, 5 -> TRACE
DEBUG = 3

##################
# Optmization
OPTIMIZE := -O2 -fopenmp

# Directories
BUILD_DIR = build

# Sources
C_SOURCES  = $(wildcard */*.c)

# Executables
COMPILE  = gcc

# Defines
C_DEFS := -DDEBUG=$(DEBUG)
ifdef SYSTEMROOT # Windows, no color prompt
C_DEFS +=
else
C_DEFS += -DLOG_USE_COLOR
endif

# Includes
C_INCLUDES := -I. -Iinc

# Compile Flags
C_FLAGS := $(C_DEFS) $(C_INCLUDES) $(OPTIMIZE) \
	-Wall -Wextra -fshort-enums -std=c11
LD_FLAGS := -lm $(OPTIMIZE)

# Objects
OBJECTS := $(addprefix $(BUILD_DIR)/,$(notdir $(C_SOURCES:.c=.o)))

vpath %.c $(sort $(dir $(C_SOURCES)))

all: $(BUILD_DIR)/$(TARGET)

$(BUILD_DIR)/%.o: %.c | $(BUILD_DIR)
	$(COMPILE) -c $< -o $@ $(C_FLAGS)
	@echo

$(BUILD_DIR)/$(TARGET): $(OBJECTS)
	$(COMPILE) -o $@ $(OBJECTS) $(LD_FLAGS)
	@echo

clean:
	rm -rf $(BUILD_DIR)

run1: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Teste 1"
	$(BUILD_DIR)/$(TARGET) 1

run2: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Teste 2"
	$(BUILD_DIR)/$(TARGET) 2

run3: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Teste 3, n=20"
	$(BUILD_DIR)/$(TARGET) 3 20
	@echo "Rodando Teste 3, n=40"
	$(BUILD_DIR)/$(TARGET) 3 40
	@echo "Rodando Teste 3, n=80"
	$(BUILD_DIR)/$(TARGET) 3 80

run41: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Rede 1"
	$(BUILD_DIR)/$(TARGET) 4 1

run42: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Rede 2"
	$(BUILD_DIR)/$(TARGET) 4 2

run43: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Rede 3"
	@echo "Isso pode levar alguns minutos"
	$(BUILD_DIR)/$(TARGET) 4 3

run44: $(BUILD_DIR)/$(TARGET)
	@echo "Rodando Rede 4"
	@echo "Isso pode levar algumas horas"
	$(BUILD_DIR)/$(TARGET) 4 4

$(BUILD_DIR):
	@echo "Creating build directory"
	@mkdir -p $@

.PHONY: all clean run1 run2 run3 run41 run42 run43 run44
