cd ../makefiles/
make clean
make config=production nthreads=8 bc=mixed precision=single particles=false dlm=true hdf=true movingframe=true bs=16 bsx=8 -j
mv simulation simulation_daint
cd ../launch/

sbatch submit_daint_fallingSamaraFixedSP
