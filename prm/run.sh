mpirun -n 1 ./maincuda 32 2048 0.01 1.0 $1 256 1
mpirun -n 2 ./maincuda 32 2048 0.01 1.0 $1 256 1
mpirun -n 4 ./maincuda 32 2048 0.01 1.0 $1 256 1
mpirun -n 1 ./maincuda 32 2048 0.01 1.0 $1 256 4
mpirun -n 2 ./maincuda 32 2048 0.01 1.0 $1 256 4
mpirun -n 4 ./maincuda 32 2048 0.01 1.0 $1 256 4
mpirun -n 1 ./maincuda 32 2048 0.01 1.0 $1 256 5
mpirun -n 2 ./maincuda 32 2048 0.01 1.0 $1 256 5
mpirun -n 4 ./maincuda 32 2048 0.01 1.0 $1 256 5

