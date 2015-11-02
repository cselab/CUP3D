#cd ../makefiles/
#make clean
#make bc=periodic config=production multiphase=false precision=double poisson=hypre -j
#cd -
#export OMP_NUM_THREADS=48
#mpirun -np 1 ../makefiles/test -minBPD 1 -maxBPD 8 -test travelingwave

cd ../makefiles/
make clean
make bc=periodic config=production multiphase=false precision=double poisson=fftw rk2=true -j
cd -
export OMP_NUM_THREADS=48
mpirun -np 1 ../makefiles/test -minBPD 1 -maxBPD 8 -minDT 1 -maxDT 1 -test travelingwave

