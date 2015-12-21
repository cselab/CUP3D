cd ../makefiles
make clean
make CC=mpic++ config=production bc=periodic precision=double bs=32 nthreads=48 rk2=true -j
cd ../launch/
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI411 -J testMPI411 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 4 -nprocsy 1 -nprocsz 1
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI141 -J testMPI141 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 1 -nprocsy 4 -nprocsz 1
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI114 -J testMPI114 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 1 -nprocsy 1 -nprocsz 4
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI221 -J testMPI221 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 2 -nprocsy 2 -nprocsz 1
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI122 -J testMPI122 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 1 -nprocsy 2 -nprocsz 2
export OMP_NUM_THREADS=48;bsub -n 192 -R 'span[ptile=48]' -W 10:00 -o testMPI212 -J testMPI212 mpirun -np 4 -pernode ../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-9 -test mpi -nprocsx 2 -nprocsy 1 -nprocsz 2