cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false -j
cd ../launch

for BPD in 32 64
#16
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Mordant_bpd${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o Mordant_${BPD} -J Mordant_${BPD} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Mordant_bpd${BPD}/FallingCylinder_Mordant_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .01 -tend 2. -rhoS 2.565 -ypos .9 -nu 0.000228609536213148 -sim falling -tdump 0.01 -lambda 1e4
#-restart -serialization /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Mordant_bpd${BPD}/FallingCylinder_Mordant_bpd${BPD}
done