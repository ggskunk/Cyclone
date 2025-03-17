# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -m64 -std=c++17 -Ofast -mssse3 -Wall -Wextra \
           -Wno-write-strings -Wno-unused-variable -Wno-deprecated-copy \
           -Wno-unused-parameter -Wno-sign-compare -Wno-strict-aliasing \
           -Wno-unused-but-set-variable \
           -march=native -mtune=native \
           -funroll-loops -ftree-vectorize -fstrict-aliasing -fno-semantic-interposition \
           -fvect-cost-model=unlimited -fno-trapping-math -fipa-ra -flto \
           -fassociative-math -fopenmp -mavx2 -mbmi2 -madx \

# Source files
SRCS = Cyclone.cpp SECP256K1.cpp Int.cpp Timer.cpp IntGroup.cpp IntMod.cpp \
       Point.cpp ripemd160_avx2.cpp sha256_avx2.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Target executable
TARGET = Cyclone

# Default target
all: fix_rdtsc $(TARGET)

# Target to replace __rdtsc with my_rdtsc
fix_rdtsc:
	find . -type f -name '*.cpp' -exec sed -i 's/__rdtsc/my_rdtsc/g' {} +

# Link the object files to create the executable and then delete .o files
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
	rm -f $(OBJS) && chmod +x $(TARGET)

# Compile each source file into an object file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	echo "Cleaning..."
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean fix_rdtsc