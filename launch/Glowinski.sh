module load gcc
cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed -j
cd ../launch

for BPD in 32
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Glowinski_bpdx${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Glowinski_bpdx${BPD}/Glowinski_bpdx${BPD} -CFL 0.1 -bpdx ${BPD} -bpdy $[${BPD}*3] -tend 10. -rhoS 1.5 -ypos .66 -nu 0.0000316227766016838 -sim falling -shape disk -radius 0.00625
done
