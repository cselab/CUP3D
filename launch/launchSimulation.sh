module load gcc
#cd ../makefiles/;make clean;make CC=mpic++ config=production poisson=multigrid -j;cd ../launch/
export OMP_NUM_THREADS=48
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMultigrid_rho.5_bpd64 -bpdx 64 -bpdy 64 -radius .025 -tend 1. -rhoS .5 -ypos .15 -nu 1.004e-3 -sim falling
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMultigrid_rho2.56_bpd64 -bpdx 64 -bpdy 64 -radius .025 -tend 1. -rhoS 2.56 -ypos .85 -nu 1.004e-3 -sim falling
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastCylinderMultigrid_rho10_bpd64 -bpdx 64 -bpdy 64 -radius .025 -tend 1. -rhoS 10 -ypos .85 -nu 1.004e-3 -sim falling

cd ../makefiles/;make clean;make config=production bc=periodic -j;cd ../launch/

bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastMovingCylinder_noMollification -bpdx 16 -bpdy 16 -radius .025 -uinf .01 -Re 100 -tend 50 -rhoS 1.00 -sim moving
#for B in 32 64 128 256
#do
#	for R in 40 100 1000
#	do
#		FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe${R}_$[${B}*32]
#		mkdir ${FOLDER_FIXED}
#		cp ../makefiles/simulation ${FOLDER_FIXED}
#		bsub -n 48 -W 168:00 ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe${R}_$[${B}*32] -bpdx ${B} -bpdy ${B} -radius .025 -uinf .01 -Re ${R} -tend 50 -rhoS 1.00 -sim fixed


#		FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe${R}_$[${B}*32]
#		mkdir ${FOLDER_MOVING}
#		cp ../makefiles/simulation ${FOLDER_MOVING}
#		bsub -n 48 -W 168:00 ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe${R}_$[${B}*32] -bpdx ${B} -bpdy ${B} -radius .025 -uinf .01 -Re ${R} -tend 50 -rhoS 1.00 -sim moving
#	done
#done

#cd ../makefiles/;make clean;make CC=mpic++ config=production poisson=split -j;cd ../launch/;
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastFallingCylinder_Split_bpd64_Mordant2560 ${SETTINGS}
#cd ../makefiles/;make clean;make CC=mpic++ config=production poisson=jacobi -j;cd ../launch/;
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastFallingCylinder_jacobi_bpd64_Mordant2560 ${SETTINGS}
#cd ../makefiles/;make clean;make CC=mpic++ config=production poisson=multigrid -j;cd ../launch/;
#bsub -n 48 -W 168:00 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FlowPastFallingCylinder_multigrid_bpd64_Mordant2560 ${SETTINGS}
