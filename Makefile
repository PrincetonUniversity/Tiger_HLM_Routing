# ==== Compiler and Flags ====
CXX := g++
CXXFLAGS := -O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG -std=c++17
#CXXFLAGS := -O0 -g -std=c++17 -fno-omit-frame-pointer -fopenmp  # For debugging

# ==== Linker Flags ====
LDFLAGS := -lboost_system -L${NETCDF_PATH}/lib64 -lnetcdf

# ==== Source Files ====
SRC := src/main.cpp \
       src/build_info.cpp \
       src/omp_info.cpp \
       src/model_setup.cpp \
       src/routing.cpp \
       src/end_info.cpp \
       src/I_O/node_info.cpp \
       src/I_O/output_series.cpp \
       src/I_O/inputs.cpp \
       src/I_O/config_loader.cpp \
       src/utils/time.cpp

# ==== Build and Binary Directories ====
BUILD_DIR := build
BIN_DIR := bin

# ==== Object files in build dir ====
OBJ := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# ==== Executable ====
BIN := $(BIN_DIR)/routing

# ==== Default Target ====
all: $(BIN)

# ==== Link executable ====
$(BIN): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ==== Compile source files ====
$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Clean ====
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean