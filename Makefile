CC          = g++
EXEC        = run
BUILD_DIR   = build

INCL        = -I./include
SRC         = src/main.cxx
OBJ         = $(BUILD_DIR)/main.o

OPTIMIZE_CC  = -fopenmp -O3 -g -Wall -Wno-unknown-pragmas -march=haswell
LINKER       = $(CC)
LINK_FLAGS   = -fopenmp

CFLAGS      = $(OPTIMIZE_CC) $(INCL)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(LINKER) $(OBJ) -o $(EXEC) $(LINK_FLAGS)

$(BUILD_DIR)/%.o: src/%.cxx | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(EXEC)

.PHONY: all clean
