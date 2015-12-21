#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 rk2=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 0 -minDT 1e-6 -maxDT 1e-6
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-6

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 rk2=true -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 0 -minDT 1e-6 -maxDT 1e-6
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-6

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 rk2=false particles=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-8
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-6

cd ../makefiles
make clean
make CC=mpic++ config=production nthreads=24 bc=periodic precision=double bs=32 rk2=true particles=false -j
cd ../launch/
export OMP_NUM_THREADS=24;mpirun -np 1 ../makefiles/test -minBPD 2 -maxBPD 4 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-8
export OMP_NUM_THREADS=24;mpirun -np 2 ../makefiles/test -minBPD 2 -maxBPD 4 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-8
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 0 -minDT 1e-8 -maxDT 1e-6



#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 1e-6 -maxDT 1e-6
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 1 -minDT 1e-9 -maxDT 1e-7

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=true -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 1e-6 -maxDT 1e-6
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -test advection -ic 1 -minDT 1e-9 -maxDT 1e-7

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=false particles=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 1e-8 -maxDT 1e-8
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 1 -minDT 1e-9 -maxDT 1e-7

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=vortex precision=double bs=32 rk2=true particles=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 32 -test advection -ic 1 -minDT 1e-6 -maxDT 1e-6
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -test advection -ic 1 -minDT 1e-9 -maxDT 1e-7