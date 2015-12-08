#cd ../makefiles
#make clean
#make config=production bc=periodic precision=double bs=32 rk2=false -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 16 -minDT 1e-9 -maxDT 1e-9 -test diffusion
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-9 -maxDT 1e-7 -test diffusion

cd ../makefiles
make clean
make config=production bc=periodic precision=double bs=32 rk2=true -j
cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 16 -minDT 1e-9 -maxDT 1e-9 -test diffusion
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 4 -minDT 1e-7 -maxDT 1e-5 -test diffusion