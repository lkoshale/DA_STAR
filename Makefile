#if clock skew detected
#find . -exec touch {} \;

all : cpu gpu

gpu : main.cpp
	nvcc -x cu main.cpp -o gm -std=c++11

cpu: main.cpp
	g++ -o m main.cpp -std=c++11

clean: gm m 
	rm -f gm m