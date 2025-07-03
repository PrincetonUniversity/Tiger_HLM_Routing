# ==== Compiler and Flags ====
CXX := g++
CXXFLAGS := -O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG -std=c++17
# CXXFLAGS := -O0 -g -std=c++17 -fno-omit-frame-pointer -fopenmp

# ==== Linker Flags ====
LDFLAGS := -lboost_system -L${NETCDF_PATH}/lib64 -lnetcdf

# ==== Project Files ====
SRC := main.cpp \
       build_info.cpp \
       omp_info.cpp \
       model_setup.cpp \
       end_info.cpp \
       I_O/node_info.cpp \
       I_O/output_series.cpp \
       I_O/inputs.cpp \
       I_O/config_loader.cpp \
       utils/time.cpp


# ==== Build directory ====
BUILD_DIR := build

# ==== Object files in build dir ====
OBJ := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRC))

# ==== Binary in build dir ====
BIN := $(BUILD_DIR)/routing

# ==== Default Target ====
all: $(BIN)

# ==== Build Binary ====
$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ==== Compile Source Files ====
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Clean Target ====
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean