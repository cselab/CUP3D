cd ../makefiles
make clean
make CC=mpic++ config=production poisson=split-fftw bc=mixed precision=double particles=false rk2=false dlm=true nthreads=12 -j
cd ../launch

for BPD in 1 2 4 8
#4 8 16 32
do
for L in .1 .5 .75 1
do
		NAME=Mordant_181215_MPI4_DLM${L}_bpd$[4*${BPD}]
		FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
		mkdir ${FOLDER}
		cp ../makefiles/simulation ${FOLDER}
		cp Mordant.sh ${FOLDER}
		export OMP_NUM_THREADS=12
		#${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -bpdz ${BPD} -radius 0.02 -tend 20. -rhoS 2.565 -ypos .9 -nu 0.000646605413200915 -sim falling -lambda ${L}
		bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} mpirun -np 4 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.01 -nprocsx 4 -nprocsy 1 -nprocsz 1 -bpdx ${BPD} -bpdy $[4*${BPD}] -bpdz $[4*${BPD}] -radius 0.02 -tend 5. -rhoS 2.56 -ypos .9 -nu 0.000646605413200915 -sim falling -dlm ${L}
	done
done
