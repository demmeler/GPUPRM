ARCHITECTURE=sm_35
#SRC=$(PWD)/src
SRC=src
BIN=bin

#CUDA_PATH=/usr/lib/nvidia-cuda-toolkit
CUDA_PATH=/usr/local/cuda-7.5
CUDA_PATH2=/share/apps/cuda-7.0
CUDA_INCLUDE= -L$(CUDA_PATH)/lib64 -L$(CUDA_PATH2)/lib64 -lcudart

NVCC_FLAGS=-lineinfo
#NVCC_FLAGS=-lineinfo -DNEW_COLLISION_KERNEL

#WRN_OPT=-Wall
WRN_OPT=

OPT=-O3
#OPT=-g -O3
#OPT=-g

#PREP=scorep
#PREP=scalasca -instrument
PREP=

#COMPILER=icpc -g
COMPILER=g++

MPICXX=mpicxx
#MPICXX=$(COMPILER) -I/usr/mpi/gcc/openmpi-1.8.4/include/openmpi/opal/mca/hwloc/hwloc191/hwloc/include -I/usr/mpi/gcc/openmpi-1.8.4/include/openmpi/opal/mca/event/libevent2021/libevent -I/usr/mpi/gcc/openmpi-1.8.4/include/openmpi/opal/mca/event/libevent2021/libevent/include -I/usr/mpi/gcc/openmpi-1.8.4/include -I/usr/mpi/gcc/openmpi-1.8.4/include/openmpi -pthread -Wl,-rpath -Wl,/usr/mpi/gcc/openmpi-1.8.4/lib64 -Wl,--enable-new-dtags -L/usr/mpi/gcc/openmpi-1.8.4/lib64 -lmpi_cxx -lmpi



all: robo arrayconfigspace $(BIN)/main $(BIN)/maincuda $(BIN)/maincuda0 $(BIN)/bench $(BIN)/benchcuda $(BIN)/benchcuda0 $(BIN)/array
	

clean:
	rm -f bin/*

arrayo:
	$(PREP) $(MPICXX) $(WRN_OPT) $(OPT) -I $(SRC)/lib -c -std=c++0x -o $(BIN)/array.o $(SRC)/array.cpp

$(BIN)/array: $(BIN)/arrayconfigspace.o $(BIN)/prmsolver.o $(BIN)/array.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)



maino:
	$(PREP) $(MPICXX) $(WRN_OPT) $(OPT) -I $(SRC)/lib -c -std=c++0x -o $(BIN)/main.o $(SRC)/main.cpp

$(BIN)/maincuda: $(BIN)/robotconfigspacecuda.o $(BIN)/prmsolver.o $(BIN)/main.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)

$(BIN)/maincuda0: $(BIN)/robotconfigspacecuda0.o $(BIN)/prmsolver.o $(BIN)/main.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)

$(BIN)/main: $(BIN)/robotconfigspace.o $(BIN)/prmsolver.o $(BIN)/main.o
	$(PREP) $(MPICXX) -o $@ $^



bencho:
	$(PREP) $(MPICXX) $(WRN_OPT) $(OPT) -I $(SRC)/lib -c -std=c++0x -o $(BIN)/bench.o $(SRC)/bench.cpp

$(BIN)/benchcuda: $(BIN)/robotconfigspacecuda.o $(BIN)/prmsolver.o $(BIN)/bench.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)

$(BIN)/benchcuda0: $(BIN)/robotconfigspacecuda0.o $(BIN)/prmsolver.o $(BIN)/bench.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)

$(BIN)/bench: $(BIN)/robotconfigspace.o $(BIN)/prmsolver.o $(BIN)/bench.o
	$(PREP) $(MPICXX) -o $@ $^ $(CUDA_INCLUDE)



# solver

prmsolver:
	$(PREP) $(MPICXX) $(WRN_OPT) $(OPT) -std=c++0x -I $(SRC)/lib -c -o $(BIN)/prmsolver.o $(SRC)/lib/prm/prmsolver.cpp


# configspaces

arrayconfigspace:
	$(PREP) g++ $(WRN_OPT) -std=c++0x -c -I $(SRC)/lib -o $(BIN)/arrayconfigspace.o $(SRC)/lib/arrayspace/arrayconfigspace.cpp


robo: robotconfigspace robotconfigspacecuda robotconfigspacecuda0
	

robotconfigspace:
	$(PREP) g++ $(WRN_OPT) $(OPT) -I $(SRC)/lib -c -o $(BIN)/robotconfigspace.o $(SRC)/lib/robotspace/robotconfigspace.cpp

robotconfigspacecuda: $(SRC)/lib/robotspace/robotconfigspace.cpp
	$(PREP) nvcc $(NVCC_FLAGS) $(OPT) -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o $(BIN)/robotconfigspacecuda.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D CUDA_CODE -D GPU_VERSION

robotconfigspacecuda0: $(SRC)/lib/robotspace/robotconfigspace.cpp
	$(PREP) nvcc $(NVCC_FLAGS) $(OPT) -x=cu -arch=$(ARCHITECTURE) -I $(SRC)/lib -c -o $(BIN)/robotconfigspacecuda0.o $(SRC)/lib/robotspace/robotconfigspace.cpp -D CUDA_CODE


# dependencies (comment out for manual compiling...)

$(BIN)/robotconfigspacecuda.o: robotconfigspacecuda
	
$(BIN)/robotconfigspacecuda0.o: robotconfigspacecuda0
	
$(BIN)/robotconfigspace.o: robotconfigspace
	
$(BIN)/prmsolver.o: prmsolver
	
$(BIN)/array.o: arrayo
	
$(BIN)/main.o: maino
	
$(BIN)/bench.o: bencho
	
$(BIN)/arrayconfigspace.o:arrayconfigspace
	


