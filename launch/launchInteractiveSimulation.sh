module load gcc

#export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastMovingCylinderRe40_interactive -bpdx 16 -bpdy 16 -radius .05 -ubody .01 -Re 40 -tend 50 -rhoS 1.01 -sim moving

#export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMultigrid_operatorReordering_Mordant1_bpd32 -bpdx 32 -bpdy 32 -radius .025 -tend 15 -rhoS 2.56 -nu 1.004e-3 -sim falling

#for B in 16 32 64
#for B in 64
#do
#	for R in 0.5 2.56 300
#for R in 2.56
#	do
#		export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderSplit_bpd${B}_rhos${R} -bpdx ${B} -bpdy ${B} -radius .025 -tend 1. -rhoS ${R} -nu 1.004e-3 -sim falling
#	done
#done

cd ../makefiles;make clean;make config=production poisson=multigrid bc=mixed -j;cd ../launch
export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_bpd64 -bpdx 64 -bpdy 64 -radius .0125 -tend 2. -rhoS 2 -ypos .85 -nu 0.0001551888 -sim falling
#export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMGSplit_rho2.56_bpd32 -bpdx 32 -bpdy 32 -radius .025 -tend 1. -rhoS 2.56 -ypos .85 -nu 1.004e-3 -sim falling

#for R in 0.01 0.1 0.5 0.9
#do
#	export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMGSplit_rho${R}_bpd16 -bpdx 16 -bpdy 16 -radius .025 -tend 1. -rhoS ${R} -ypos 0.15 -nu 1.004e-3 -sim falling
#done

#for R in 1.01 2.56 10 20 50 100 300
#do
#	export OMP_NUM_THREADS=48;../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMGSplit_rho${R}_bpd16 -bpdx 16 -bpdy 16 -radius .025 -tend 1. -rhoS ${R} -ypos 0.85 -nu 1.004e-3 -sim falling
#done