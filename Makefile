#if clock skew detected
#find . -exec touch {} \;
# Makefile for Writing Make Files Example

# *****************************************************
# Variables to control Makefile operation

NVCC = nvcc
NVCCFLAGS = -x cu 
CXX = g++
CXXFLAGS = -Wall 

# ****************************************************
# Targets needed to bring the executable up to date

all: gpu cpu

gpu: main.cpp
	$(NVCC) $(NVCCFLAGS) main.cpp -o gpu -std=c++11

cpu: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o cpu -std=c++11
	
clean: gpu cpu
	rm -f gpu cpu
