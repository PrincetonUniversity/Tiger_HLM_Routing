# ==== Compiler and Flags ====
CXX := g++
CXXFLAGS := -O3 -march=native -ffast-math -funroll-loops -fopenmp -DNDEBUG -std=c++17
# CXXFLAGS := -O0 -g -std=c++17 -fno-omit-frame-pointer -fopenmp

# ==== Linker Flags ====
LDFLAGS := -lboost_system -L${NETCDF_PATH}/lib64 -lnetcdf

# ==== Project Files ====
SRC := routing_test.cpp\
	   I_O/node_info.cpp \
	   I_O/output_series.cpp \
	   I_O/inputs.cpp
	   
OBJ := $(SRC:.cpp=.o)
BIN := routing

# ==== Default Target ====
all: $(BIN)

# ==== Build Binary ====
$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# ==== Compile Source Files ====
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


# ==== Clean Target ====
clean:
	rm -f $(OBJ) $(BIN)

# ==== Future Targets ====
# test: Run unit tests
# bench: Benchmark performance
# debug: Compile with debugging symbols
# install: Install to system
# doc: Generate documentation

.PHONY: all clean
