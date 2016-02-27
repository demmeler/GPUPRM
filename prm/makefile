ARCHITECTURE=sm_35

all:
#	g++ -o test test.cpp
#	g++ -o test2 test2.cpp
#	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test2 test2.cpp
#	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test4 test4.cpp
#	g++ -o colltest colltest.cpp
#	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp
#	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp
#	g++ -o colltest2 colltest2.cpp
#	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltest2 colltest2.cpp
#	g++ -O3 -std=c++11 -o test5 test5.cpp
	g++ -g -O3 -std=c++11 -o testvertexlist testvertexlist.cpp
#	g++ -g -std=c++11 -o testvertexlist testvertexlist.cpp
#	g++ -g -O3 -std=c++11 -o testvertexlist2 testvertexlist2.cpp
#	g++ -g -std=c++11 -o testvertexlist2 testvertexlist2.cpp

test6: arrayconfigspace
	g++ -g -c -std=c++11 -o test6.o test6.cpp
	g++ -g -o test6 test6.o arrayconfigspace.o

test3: arrayconfigspace
	g++ -c -o test3.o test3.cpp
	g++ -o test3 test3.o arrayconfigspace.o

arrayconfigspace:
	g++ -c -g -o arrayconfigspace.o lib/arrayconfigspace.cpp
