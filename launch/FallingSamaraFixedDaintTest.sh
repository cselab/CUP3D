cd ../makefiles/
make clean
make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 bsx=4 -j
mv simulation simulation_bsx4
make clean
make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 bsx=8 -j
mv simulation simulation_bsx8
make clean
make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 -j
mv simulation simulation_bsx16
cd ../launch/

sbatch submit_daint_fallingSamaraFixed_4
sbatch submit_daint_fallingSamaraFixed_8
sbatch submit_daint_fallingSamaraFixed_16