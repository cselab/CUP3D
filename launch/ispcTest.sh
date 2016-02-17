cd ../makefiles/
make clean
make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 ispc=true rk2=false -j
#make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=16 bsx=16 ispc=true rk2=false -j
cd ../launch/

export OMP_NUM_THREADS=8
aprun -n 1 -N 1 -d 8 ../makefiles/simulation -file ISPCTest -sim falling -tend 1000 -nprocsx 1 -nprocsy 1 -nprocsz 1 -bpdx 4 -bpdy 4 -bpdz 4 -CFL .01 -shape sphere -radius 0.02 -rhoS 2.56 -lambda 1 -dlm 1 -nu 0.00001 -isosurface .1 -fdump 1000
#sbatch submit_daint_fallingSamaraFixed