cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false dlm=false -j
cd ../launch

BPD=16

for L in 1e2 1e4 1e6
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Mordant_large_Lambda${L}_bpd${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
#export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o Mordant_Lambda${L}_${BPD} -J Mordant_Lambda${L}_${BPD} mpirun -np 32 ${FOLDER}/simulation -file ${FOLDER}/Mordant_Lambda${L}_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .01 -tend 2. -rhoS 2.565 -ypos .9 -nu 0.000228609536213148 -sim falling -tdump 0.01 -lambda ${L}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o Mordant_Lambda${L}_${BPD} -J Mordant_Lambda${L}_${BPD} mpirun -np 32 ${FOLDER}/simulation -file ${FOLDER}/Mordant_Lambda${L}_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .02 -tend 2. -rhoS 2.565 -ypos .9 -nu 0.000646605413200915 -sim falling -tdump 0.01 -lambda ${L}
done


#DLM
cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false dlm=true -j
cd ../launch

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Mordant_large_DLM_bpd${BPD}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o Mordant_DLM_${BPD} -J Mordant_DLM_${BPD} mpirun -np 32 ${FOLDER}/simulation -file ${FOLDER}/Mordant_DLM_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .01 -tend 2. -rhoS 2.565 -ypos .9 -nu 0.000228609536213148 -sim falling -tdump 0.01