#cd ../makefiles
#make clean
#make config=production bc=periodic precision=single bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -solver 0 -ic 0 -minDT 1 -maxDT 1 > results_pressure_FFTW_stencil_Poisson
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -solver 0 -ic 1 -minDT 1 -maxDT 1 > results_pressure_FFTW_stencil_Vel
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -solver 1 -ic 0 -minDT 1 -maxDT 1 > results_pressure_FFTW_spectral_Poisson
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -solver 1 -ic 1 -minDT 1 -maxDT 1 > results_pressure_FFTW_spectral_Vel

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 2 -ic 0 -minDT 1 -maxDT 1 > results_pressure_MGconst_Poisson
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 2 -ic 1 -minDT 1 -maxDT 1 > results_pressure_MGconst_Vel

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 0 -minDT 1 -maxDT 1 > results_pressure_MG_Poisson
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 1 -minDT 1 -maxDT 1 > results_pressure_MG_Vel

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=single bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;time ../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 2 -ic 2 -minDT 1 -maxDT 1
#export OMP_NUM_THREADS=48;time ../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 0 -ic 2 -minDT 1 -maxDT 1


cd ../makefiles
make clean
make config=debug poisson=fftw bc=mixed precision=double bs=32 -j
cd ../launch/
#export OMP_NUM_THREADS=48;time ../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 2 -ic 2 -minDT 1 -maxDT 1
export OMP_NUM_THREADS=48;time ../makefiles/test -minBPD 2 -maxBPD 16 -test pressure -solver 0 -ic 2 -minDT 1 -maxDT 1