# ==== Compiler and Flags ====
CXX := mpiicpx
NVCC := nvcc
CXXFLAGS := -O3 -ipo -fp-model fast=2 -qopenmp -fma \
            -xSapphireRapids -qopt-report-phase=vec -qopt-prefetch \
            -Rpass=loop-vectorize -Rpass=inline -DNDEBUG -std=c++17

# ==== Build and Binary Directories ====
BUILD_DIR := build
BIN_DIR := bin

# GPU support: build with `make USE_GPU=1` to enable level-0 GPU solving
# Without USE_GPU=1, the binary is identical to the original CPU-only build.
USE_GPU ?= 0
ifeq ($(USE_GPU),1)
CXXFLAGS += -DUSE_GPU_LEVEL0
NVCCFLAGS := -O3 -std=c++17 --expt-relaxed-constexpr \
                 -Xcompiler -fopenmp \
                 -I$(BOOST_ROOT)/include -I$(NETCDF_PATH)/include
LDFLAGS_GPU := -L$(CUDA_HOME)/lib64 -lcudart
GPU_OBJ := $(BUILD_DIR)/models/level0_gpu.o
else
LDFLAGS_GPU :=
GPU_OBJ :=
endif

# ==== Linker Flags ====
LDFLAGS := -lboost_system -L${NETCDF_PATH}/lib64 -lnetcdf $(LDFLAGS_GPU)

# ==== Source Files ====
SRC := src/main.cpp \
       src/build_info.cpp \
       src/omp_info.cpp \
       src/model_setup.cpp \
       src/dependency_graph.cpp \
       src/partition.cpp \
       src/boundary_exchange.cpp \
       src/routing.cpp \
       src/end_info.cpp \
       src/I_O/node_info.cpp \
       src/I_O/output_series.cpp \
       src/I_O/inputs.cpp \
       src/I_O/config_loader.cpp \
       src/utils/time.cpp \
       src/utils/level_timing.cpp

# ==== Object files in build dir ====
OBJ := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRC))
OBJ += $(GPU_OBJ)

# ==== Executable ====
BIN := $(BIN_DIR)/routing

# ==== Default Target ====
all: $(BIN)
	@echo "Build successful: $(BIN)"

# ==== Link executable ====
$(BIN): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ==== Compile C++ source files ====
$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ==== Compile CUDA source file (only when USE_GPU=1) ====
ifeq ($(USE_GPU),1)
$(BUILD_DIR)/models/level0_gpu.o: src/models/level0_gpu.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
endif

# ==== Clean ====
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean