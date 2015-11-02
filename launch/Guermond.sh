cd ../makefiles/
make clean
make poisson=hypre bc=mixed config=production multiphase=true particles=false -j
cd -
export OMP_NUM_THREADS=48
#mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/rti_Guermond_veryHR_Re1000 -CFL 0.05 -bpdx 64 -bpdy 256 -tend 10. -sim rti -tdump 0.02
#bsub -W 168:00 -n 48 -o Guermond_HR -J Guermond_HR mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/rti_Guermond_HR_Re1000 -CFL 0.05 -bpdx 32 -bpdy 128 -tend 10. -sim rti -tdump 0.01
bsub -W 168:00 -n 48 -o Guermond_FD_MR -J Guermond_FD_MR mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/rti_Guermond_FD_MR_Re1000 -CFL 0.05 -bpdx 16 -bpdy 64 -tend 10. -sim rti -tdump 0.01
