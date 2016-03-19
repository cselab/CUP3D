cd ../makefiles/
make clean
make config=production nthreads=12 bc=mixed precision=double particles=false rk2=false dlm=true hdf=true movingframe=true bs=16 bsx=8 sz=0.963 -j
mv simulation simulation_dora_dp
cd ../launch/

sbatch submit_dora_fallingSamaraFixed
