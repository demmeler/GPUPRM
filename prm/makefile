ARCHITECTURE=sm_20
SRC=src
#CUDA_PATH=/usr/lib/nvidia-cuda-toolkit
CUDA_PATH=/usr/local/cuda-7.5
CUDA_PATH2=/share/apps/cuda-5.5
CUDA_INCLUDE= -L$(CUDA_PATH)/lib64 -L$(CUDA_PATH2)/lib64 -lcudart


all: test10 test10cuda


test6: arrayconfigspace
	g++ -g -c -std=c++11 -o test6.o $(SRC)/test6.cpp
	g++ -g -o test6 test6.o arrayconfigspace.o

test9o:
	g++ -O3 -g -c -std=c++0x -o test9.o $(SRC)/test9.cpp

test9cuda: robotconfigspacecuda.o test9o
	g++ -g -o test9cuda test9.o robotconfigspacecuda.o $(CUDA_INCLUDE)

test9: robotconfigspacecuda0.o test9o
	g++ -g -o test9 test9.o robotconfigspacecuda0.o $(CUDA_INCLUDE)


test10o:
	mpicxx -O3 -g -c -std=c++0x -o test10.o $(SRC)/test10.cpp

test10cuda: robotconfigspacecuda.o test10o
	mpicxx -g -o test10cuda test10.o robotconfigspacecuda.o $(CUDA_INCLUDE)

test10: robotconfigspacecuda0.o test10o
	mpicxx -g -o test10 test10.o robotconfigspacecuda0.o $(CUDA_INCLUDE)



# configspaces


arrayconfigspace:
	g++ -std=c++11 -c -g -o arrayconfigspace.o $(SRC)/lib/arrayconfigspace.cpp


robo: robotconfigspacecuda robotconfigspacecuda0

#robotconfigspace:
#	g++ -std=c++11 -c -g -o robotconfigspace.o $(SRC)/lib/robotconfigspace.cpp

robotconfigspacecuda: $(SRC)/lib/robotconfigspace.cpp
	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o robotconfigspacecuda.o $(SRC)/lib/robotconfigspace.cpp -D CUDA_IMPLEMENTATION

robotconfigspacecuda0: $(SRC)/lib/robotconfigspace.cpp
	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -c -o robotconfigspacecuda0.o $(SRC)/lib/robotconfigspace.cpp

