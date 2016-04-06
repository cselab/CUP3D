cd ../makefiles/
make clean
make config=production nthreads=8 bc=mixed precision=single particles=false rk2=true dlm=true hdf=true movingframe=true j0=false bs=16 bsx=8 -j
mv simulation simulation_daint
cd ../launch/

sbatch submit_daint_fallingSamaraFixedSP
