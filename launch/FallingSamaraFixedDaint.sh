cd ../makefiles/
make clean
#make config=production nthreads=2 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bsx=8 -j
make config=production nthreads=8 bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true bs=32 bsx=8 -j
cd ../launch/

#aprun -n 4 -N 4 -d 2 ../makefiles/simulation -file test -sim falling -tend 1000 -nprocsx 4 -nprocsy 1 -nprocsz 1 -bpdx 4 -bpdy 4 -bpdz 4 -CFL .1 -shape samara -rhoS 250. -lambda 1 -dlm 1 -nu 0.00001 -isosurface .01 -tdump .01 -fdump 500
sbatch submit_daint_fallingSamaraFixed