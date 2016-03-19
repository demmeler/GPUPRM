ARCHITECTURE=sm_35
SRC=src
CUDA_PATH=/usr/lib/nvidia-cuda-toolkit
CUDA_INCLUDE= -L$(CUDA_PATH)/lib64 -lcudart

test:
	g++ -o test $(SRC)/test.cpp

test2:
#	g++ -o test2 $(SRC)/test2.cpp
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test2 $(SRC)/test2.cpp

test4:
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o test4 $(SRC)/test4.cpp

colltest:
	g++ -o colltest $(SRC)/colltest.cpp

colltestcuda:
	nvcc -O3 -x=cu -arch=$(ARCHITECTURE) -o colltestcuda $(SRC)/colltestcuda.cpp
#	nvcc -g -G -x=cu -arch=$(ARCHITECTURE) -o colltestcuda $(SRC)/colltestcuda.cpp

colltest2:
#	g++ -o colltest2 $(SRC)/colltest2.cpp
	nvcc -std=c++11 -g -G -x=cu -arch=$(ARCHITECTURE) -o colltest2 $(SRC)/colltest2.cpp

test5:
	g++ -O3 -std=c++11 -o test5 $(SRC)/test5.cpp

testvertexlist:
	g++ -g -O3 -std=c++11 -o testvertexlist $(SRC)/testvertexlist.cpp
#	g++ -g -std=c++11 -o testvertexlist $(SRC)/testvertexlist.cpp
testvertexlist2:
	g++ -g -O3 -std=c++11 -o testvertexlist2 $(SRC)/testvertexlist2.cpp
#	g++ -g -std=c++11 -o testvertexlist2 $(SRC)/testvertexlist2.cpp

test6: arrayconfigspace
	g++ -g -c -std=c++11 -o test6.o $(SRC)/test6.cpp
	g++ -g -o test6 test6.o arrayconfigspace.o

test7: robotconfigspace
	g++ -g -c -std=c++11 -o test7.o $(SRC)/test7.cpp
	g++ -g -o test7 test7.o robotconfigspace.o

test8: robotconfigspace
	g++ -g -c -std=c++11 -o test8.o $(SRC)/test8.cpp
	g++ -g -o test8 test8.o robotconfigspace.o

test9: robotconfigspace
	g++ -g -c -std=c++11 -o test9.o $(SRC)/test9.cpp
	g++ -g -o test9 test9.o robotconfigspace.o

test9cuda: robotconfigspacecuda
	g++ -O3 -g -c -std=c++11 -o test9.o $(SRC)/test9.cpp
	g++ -g -o test9cuda test9.o robotconfigspacecuda.o $(CUDA_INCLUDE)

test8cuda: robotconfigspacecuda
	nvcc -std=c++11 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o test8cuda.o $(SRC)/test8cuda.cpp
#	nvcc -std=c++11 -g -G -x=cu -arch=$(ARCHITECTURE) -o test8cuda test8cuda.o robotconfigspacecuda.o
	g++ -g -o test8cuda test8cuda.o robotconfigspacecuda.o -L"/usr/local/cuda-7.5"/lib64 -lcudart

test3: arrayconfigspace
	g++ -c -o test3.o $(SRC)/test3.cpp
	g++ -o test3 test3.o arrayconfigspace.o

arrayconfigspace:
	g++ -std=c++11 -c -g -o arrayconfigspace.o $(SRC)/lib/arrayconfigspace.cpp

robotconfigspace:
	g++ -std=c++11 -c -g -o robotconfigspace.o $(SRC)/lib/robotconfigspace.cpp

robotconfigspacecuda:
#	nvcc -std=c++11 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o robotconfigspacecuda.o $(SRC)/lib/robotconfigspace.cpp
#	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o robotconfigspacecuda.o $(SRC)/lib/robotconfigspace.cpp
	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o robotconfigspacecuda.o $(SRC)/lib/robotconfigspace.cpp -D CUDA_IMPLEMENTATION


