cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double rk2=true particles=false dlm=true -j
#make config=production poisson=split-fftw bc=periodic precision=double rk2=true particles=false dlm=true -j
cd ../launch
for C in 0.01 0.1
#0.05 0.2
do
	for L in 1e5
#1e6
#1e4
	do
		for BPD in 8 16
#8 16 32
		do
			NAME=Falling_DCT_Re40_0211_DLM10_bpd${BPD}_CFL${C}
			FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
			mkdir ${FOLDER}
			cp ../makefiles/simulation ${FOLDER}
            cp FallingRe40.sh ${FOLDER}
            cd ${FOLDER}
            export OMP_NUM_THREADS=48
            bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ./simulation -file ${FOLDER}/${NAME} -CFL ${C} -LCFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 2. -rhoS 2. -ypos .9 -nu 0.000846515308813225 -sim falling -lambda ${L} -DLM 10
            cd -
		done
	done
done


#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=single particles=false rk2=false -j
#cd ../launch
#for C in 0.05
#do
#    for L in 1e6
#    do
#        for BPD in 8 16 32
#        do
#            NAME=Falling_MG_Re40_1910_SP_Lambda${L}_bpd${BPD}_CFL${C}
#            FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
#            mkdir ${FOLDER}
#            cp ../makefiles/simulation ${FOLDER}
#            export OMP_NUM_THREADS=48
#            bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -LCFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 3. -rhoS 2. -ypos .9 -nu 0.000846515308813225 -sim falling -tdump 0.025 -lambda ${L}
#        done
#    done
#done
