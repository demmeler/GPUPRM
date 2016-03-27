ARCHITECTURE=sm_20
#SRC=$(PWD)/src
SRC=src
#CUDA_PATH=/usr/lib/nvidia-cuda-toolkit
CUDA_PATH=/usr/local/cuda-7.5
CUDA_PATH2=/share/apps/cuda-5.5
CUDA_INCLUDE= -L$(CUDA_PATH)/lib64 -L$(CUDA_PATH2)/lib64 -lcudart



all: robo arrayconfigspace main maincuda test10 test10cuda

test10o:
	mpicxx -O3 -g -I $(SRC)/lib -c -std=c++0x -o test10.o $(SRC)/test10.cpp

test10cuda: robotconfigspacecuda.o prmsolver.o test10.o
	mpicxx -g -o test10cuda test10.o prmsolver.o robotconfigspacecuda.o $(CUDA_INCLUDE)

test10: robotconfigspacecuda0.o prmsolver.o test10.o
	mpicxx -g -o test10 test10.o prmsolver.o robotconfigspacecuda0.o $(CUDA_INCLUDE)


maino:
	mpicxx -O3 -g -I $(SRC)/lib -c -std=c++0x -o main.o $(SRC)/main.cpp

maincuda: robotconfigspacecuda.o prmsolver.o main.o
	mpicxx -g -o maincuda main.o prmsolver.o robotconfigspacecuda.o $(CUDA_INCLUDE)

main: robotconfigspacecuda0.o prmsolver.o main.o
	mpicxx -g -o main main.o prmsolver.o robotconfigspacecuda0.o $(CUDA_INCLUDE)



# solver

prmsolver:
	mpicxx -O3 -g -std=c++0x -I $(SRC)/lib -c -o prmsolver.o $(SRC)/lib/prm/prmsolver.cpp


# configspaces

arrayconfigspace:
	g++ -std=c++11 -c -g -I $(SRC)/lib -o arrayconfigspace.o $(SRC)/lib/arrayspace/arrayconfigspace.cpp


robo: robotconfigspacecuda robotconfigspacecuda0

robotconfigspace:
	g++ -O3 -g -I $(SRC)/lib -c -o robotconfigspace.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D NO_CUDA

robotconfigspacecuda: $(SRC)/lib/robotspace/robotconfigspace.cpp
	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o robotconfigspacecuda.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D CUDA_IMPLEMENTATION

robotconfigspacecuda0: $(SRC)/lib/robotspace/robotconfigspace.cpp
	nvcc -O3 -g -G -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o robotconfigspacecuda0.o $(SRC)/lib/robotspace/robotconfigspace.cpp


# dependancies (comment out for manual compiling...)

robotconfigspacecuda.o: robotconfigspacecuda
robotconfigspacecuda0.o: robotconfigspacecuda0
prmsolver.o: prmsolver
test10.o: test10o
main.o: maino
