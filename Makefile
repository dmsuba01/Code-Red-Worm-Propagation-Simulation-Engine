# ==============================================================================
# Makefile for Code Red Worm Propagation Simulation Engine (OpenMP)
# Authors: Bhavashruthi M & D M Subashri (Team T05 / Project P20)
# ==============================================================================

CXX = g++
CXXFLAGS = -O3 -Wall -std=c++17 -fopenmp
LDFLAGS = -fopenmp

TARGET = worm_simulation
SRC = worm_simulation.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET) 100

clean:
	rm -f $(TARGET) $(TARGET).exe *.o

