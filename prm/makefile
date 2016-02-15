ARCHITECTURE=sm_11

all:
#	g++ -o test test.cpp
#	g++ -o test2 test2.cpp
#	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test2 test2.cpp
#	g++ -o colltest colltest.cpp
#	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp
#	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp

test3: arrayconfigspace
	g++ -c -o test3.o test3.cpp
	g++ -o test3 test3.o arrayconfigspace.o

arrayconfigspace:
	g++ -c -o arrayconfigspace.o lib/arrayconfigspace.cpp
