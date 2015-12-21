cd ../makefiles
make clean
make CC=mpic++ config=production poisson=split-fftw bc=mixed precision=double particles=false rk2=false dlm=true movingframe=true -j
cd ../launch

for BPD in 2 4 8 16
#4 8 16 32
do
for L in .1 .5 .75 1
do
		NAME=Mordant_Fixed_151215_DLM${L}_bpd${BPD}
		FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
		mkdir ${FOLDER}
		cp ../makefiles/simulation ${FOLDER}
		cp Mordant.sh ${FOLDER}
		export OMP_NUM_THREADS=48
		#${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -bpdz ${BPD} -radius 0.02 -tend 20. -rhoS 2.565 -ypos .9 -nu 0.000646605413200915 -sim falling -lambda ${L}
		bsub -n 48 -W 40:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.01 -bpdx ${BPD} -bpdy ${BPD} -bpdz ${BPD} -radius 0.02 -tend 5. -rhoS 2.56 -ypos .9 -nu 0.000646605413200915 -sim falling -dlm ${L}
	done
done
