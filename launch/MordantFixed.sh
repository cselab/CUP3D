P=$1

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/

cd ../makefiles
make clean
make CC=mpic++ config=production poisson=split-fftw bc=mixed precision=double particles=false rk2=false dlm=true movingframe=true nthreads=48 -j
cd ../launch

for BPD in 1 2 4 8 16 32
#1 2 4 8
#4 8 16 32
do
for L in .1 .5 .75 1
do
		NAME=MordantFixed_281215_MPI${P}_DLM${L}_bpd$[4*${BPD}]
		FOLDER=${BASEPATH}${NAME}
		mkdir ${FOLDER}
		cp ../makefiles/simulation ${FOLDER}
		cp MordantFixed.sh ${FOLDER}
		cd ${FOLDER}
		export OMP_NUM_THREADS=48
		#${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -bpdz ${BPD} -radius 0.02 -tend 20. -rhoS 2.565 -ypos .9 -nu 0.000646605413200915 -sim falling -lambda ${L}
		bsub -n $[48*${P}] -R 'span[ptile=48]' -W 40:00 -o ${NAME} -J ${NAME} mpirun -np ${P} -pernode ./simulation -file ${NAME} -CFL 0.01 -nprocsx ${P} -nprocsy 1 -nprocsz 1 -bpdx ${BPD} -bpdy $[${P}*${BPD}] -bpdz $[${P}*${BPD}] -radius 0.02 -tend 5. -rhoS 2.56 -ypos .9 -nu 0.000646605413200915 -sim falling -dlm ${L}
		cd -
	done
done
