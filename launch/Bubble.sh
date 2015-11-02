cd ../makefiles/
make clean
make poisson=hypre bc=mixed config=production multiphase=true particles=true -j
cd -
export OMP_NUM_THREADS=48
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp_testP -LCFL 0.1 -bpdx 8 -bpdy 8 -sim bubble -rhoS 1000 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 10000 -split
export OMP_NUM_THREADS=48;mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/Bubble -CFL 0.1 -bpdx 16 -bpdy 16 -radius .0125 -tend 2. -rhoS 1000 -ypos .85 -nu 0.001551888 -sim bubble -fdump 100 -split
#for I in 2 4 8 16 32
#do
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS0.01 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 0.01 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS0.1 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 0.1 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS0.5 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 0.5 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS0.9 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 0.9 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS1.0 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 1.0 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS1.1 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 1.1 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS2.0 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 2.0 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS10.0 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 10.0 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp${I}_rhoS100.0 -LCFL 0.1 -bpdx ${I} -bpdy ${I} -sim bubble -rhoS 100.0 -radius 0.05 -fdump 100 -CFL 0.1 -nsteps 100
#done