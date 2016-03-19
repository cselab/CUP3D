cd ../makefiles/
make clean
make config=production nthreads=12 bc=mixed precision=single particles=false rk2=false dlm=true hdf=true movingframe=true j0=false bs=16 bsx=8 -j
mv simulation simulation_dora
cd ../launch/

sbatch submit_dora_fallingSamaraFixedSP
