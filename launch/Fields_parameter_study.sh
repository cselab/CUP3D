#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=double particles=false -j
#cd ../launch

#for BPD in 32
#do
#	for H in .0125 .00625 .003125
#	do
#		for N in 0.0000001 0.000001 0.00001 0.0001
#		do
#			for D in 1.01 1.1 1.5 2.
#			do
#				BASENAME=Fields_paramTest_bpd${BPD}_H${H}_Nu${N}_rhoS${D}
#				FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
#				mkdir ${FOLDER}
#				cp ../makefiles/simulation ${FOLDER}
#				export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME} -J ${BASENAME} mpirun -np 32 ${FOLDER}/simulation -file ${FOLDER}/${BASENAME} -CFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS ${D} -ypos .9 -nu ${N} -sim falling -shape ellipse -semiAxisX ${H} -semiAxisY .025 -angle 0.79 -tdump 0.025
#			done
#		done
#	done
#done

for BPD in 64
do
	for H in .00625 .003125
	do
		for N in 0.000001 0.00001 0.0001
		do
			for D in 1.01 1.1 1.5 2.
			do
				BASENAME=Fields_paramTest_bpd${BPD}_H${H}_Nu${N}_rhoS${D}
				FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
#mkdir ${FOLDER}
#cp ../makefiles/simulation ${FOLDER}
				export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME} -J ${BASENAME} mpirun -np 32 ${FOLDER}/simulation -file ${FOLDER}/${BASENAME} -CFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS ${D} -ypos .9 -nu ${N} -sim falling -shape ellipse -semiAxisX ${H} -semiAxisY .025 -angle 0.79 -tdump 0.05 -serialization ${FOLDER_1}/${BASENAME}${BPD}_ic1 -restart

				for N in {1..10}
				do
					export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -w "ended(${BASENAME}${BPD}_ic1_$[${N}-1])" -o ${BASENAME}${BPD}_ic1_$N -J ${BASENAME}${BPD}_ic1_$N mpirun -np 32 ${FOLDER_1}/simulation -serialization ${FOLDER_1}/${BASENAME}${BPD}_ic1 -restart -sim falling
				done
			done
		done
	done
done