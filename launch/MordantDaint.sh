cd ../makefiles/
make clean
make config=production nthreads=8 bc=mixed precision=single rk2=false particles=false dlm=true hdf=true movingframe=true j0=false bs=32 bsx=16 -j
mv simulation simulation_mordant
cd ../launch/

sbatch submit_Mordant
#sbatch submit_MordantTungsten