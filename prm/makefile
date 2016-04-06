ARCHITECTURE=sm_20
#SRC=$(PWD)/src
SRC=src
#CUDA_PATH=/usr/lib/nvidia-cuda-toolkit
CUDA_PATH=/usr/local/cuda-7.5
CUDA_PATH2=/share/apps/cuda-5.5
CUDA_INCLUDE= -L$(CUDA_PATH)/lib64 -L$(CUDA_PATH2)/lib64 -lcudart
#WRN_OPT=-Wall
WRN_OPT=
OPT=-O3
#OPT=


all: robo arrayconfigspace main maincuda maincuda0 bench benchcuda array


arrayo:
	mpicxx $(WRN_OPT) $(OPT) -g -I $(SRC)/lib -c -std=c++0x -o array.o $(SRC)/array.cpp

array: arrayconfigspace.o prmsolver.o array.o
	mpicxx -g -o array array.o prmsolver.o arrayconfigspace.o $(CUDA_INCLUDE)



maino:
	mpicxx $(WRN_OPT) $(OPT) -g -I $(SRC)/lib -c -std=c++0x -o main.o $(SRC)/main.cpp

maincuda: robotconfigspacecuda.o prmsolver.o main.o
	mpicxx -g -o maincuda main.o prmsolver.o robotconfigspacecuda.o $(CUDA_INCLUDE)

maincuda0: robotconfigspacecuda0.o prmsolver.o main.o
	mpicxx -g -o maincuda0 main.o prmsolver.o robotconfigspacecuda0.o $(CUDA_INCLUDE)

main: robotconfigspace.o prmsolver.o main.o
	mpicxx -g -o main main.o prmsolver.o robotconfigspace.o



bencho:
	mpicxx $(WRN_OPT) $(OPT) -g -I $(SRC)/lib -c -std=c++0x -o bench.o $(SRC)/bench.cpp

benchcuda: robotconfigspacecuda.o prmsolver.o bench.o
	mpicxx -g -o benchcuda bench.o prmsolver.o robotconfigspacecuda.o $(CUDA_INCLUDE)

bench: robotconfigspacecuda0.o prmsolver.o bench.o
	mpicxx -g -o bench bench.o prmsolver.o robotconfigspacecuda0.o $(CUDA_INCLUDE)




# solver

prmsolver:
	mpicxx $(WRN_OPT) $(OPT) -g -std=c++0x -I $(SRC)/lib -c -o prmsolver.o $(SRC)/lib/prm/prmsolver.cpp


# configspaces

arrayconfigspace:
	g++ $(WRN_OPT) -std=c++0x -c -g -I $(SRC)/lib -o arrayconfigspace.o $(SRC)/lib/arrayspace/arrayconfigspace.cpp


robo: robotconfigspace robotconfigspacecuda robotconfigspacecuda0

robotconfigspace:
	g++ $(WRN_OPT) $(OPT) -g -I $(SRC)/lib -c -o robotconfigspace.o $(SRC)/lib/robotspace/robotconfigspace.cpp

robotconfigspacecuda: $(SRC)/lib/robotspace/robotconfigspace.cpp
	nvcc $(OPT) -g -G -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o robotconfigspacecuda.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D CUDA_CODE -D GPU_VERSION

robotconfigspacecuda0: $(SRC)/lib/robotspace/robotconfigspace.cpp
	nvcc $(OPT) -g -G -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o robotconfigspacecuda0.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D CUDA_CODE


# dependencies (comment out for manual compiling...)

robotconfigspacecuda.o: robotconfigspacecuda
robotconfigspacecuda0.o: robotconfigspacecuda0
robotconfigspace.o: robotconfigspace
prmsolver.o: prmsolver
array.o: arrayo
main.o: maino
bench.o: bencho
arrayconfigspace.o:arrayconfigspace
