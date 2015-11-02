cd ../makefiles/
make clean
make poisson=hypre bc=periodic multiphase=true config=production -j
cd -
export OMP_NUM_THREADS=48
mpirun -np 1 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/jet_rhos${R} -CFL 0.1 -bpdx 8 -bpdy 8 -tend 10. -rhoS ${R} -sim jet -tdump 0.1
