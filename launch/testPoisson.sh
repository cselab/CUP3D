cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;mpirun -np 1 ../makefiles/test -minBPD 8 -maxBPD 64 -minDT 1 -maxDT 1 -test poisson -ic 0
export OMP_NUM_THREADS=48;mpirun -np 1 ../makefiles/test -minBPD 8 -maxBPD 64 -minDT 1 -maxDT 1 -test poisson -ic 1
