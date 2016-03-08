cd ../makefiles/
make clean
make config=production nthreads=24 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 bsx=4 -j
cd ../launch/

sbatch submit_dora_fallingSamaraFixed_restartTest
