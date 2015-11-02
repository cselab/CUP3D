cd ../makefiles
make clean
make config=production bc=mixed precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -test bc -minBPD 2 -maxBPD 2