ARCHITECTURE=sm_35

test:
	g++ -o test test.cpp

test2:
#	g++ -o test2 test2.cpp
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test2 test2.cpp

test4:
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test4 test4.cpp

colltest:
	g++ -o colltest colltest.cpp

colltestcuda:
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp
#	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltestcuda colltestcuda.cpp

colltest2:
#	g++ -o colltest2 colltest2.cpp
	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltest2 colltest2.cpp

test5:
	g++ -O3 -std=c++11 -o test5 test5.cpp

testvertexlist:
	g++ -g -O3 -std=c++11 -o testvertexlist testvertexlist.cpp
#	g++ -g -std=c++11 -o testvertexlist testvertexlist.cpp
testvertexlist2:
	g++ -g -O3 -std=c++11 -o testvertexlist2 testvertexlist2.cpp
#	g++ -g -std=c++11 -o testvertexlist2 testvertexlist2.cpp

test6: arrayconfigspace
	g++ -g -c -std=c++11 -o test6.o test6.cpp
	g++ -g -o test6 test6.o arrayconfigspace.o

test3: arrayconfigspace
	g++ -c -o test3.o test3.cpp
	g++ -o test3 test3.o arrayconfigspace.o

arrayconfigspace:
	g++ -std=c++11 -c -g -o arrayconfigspace.o lib/arrayconfigspace.cpp
