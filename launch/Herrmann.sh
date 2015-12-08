BASENAME=RTI3D_251115
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
mkdir ${FOLDER}

cd ../makefiles/
make clean
#does not seem to work with RK2
make poisson=split-fftw bc=mixed config=production multiphase=true particles=false rk2=false -j
cd -
for BPD in 8 16 32
#2 4
do
    cp ../makefiles/simulation ${FOLDER}/simulation_BPD${BPD}_FD
	export OMP_NUM_THREADS=48
	bsub -n 48 -W 168:00 -R "rusage[mem=4096]" -o ${BASENAME}_BPD${BPD}_FD -J ${BASENAME}_BPD${BPD}_FD ${FOLDER}/simulation_BPD${BPD}_FD -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/Herrmann_FD_${BPD} -CFL 0.5 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -bpdz ${BPD} -tend 1.25 -sim rti -fdump $[20*${BPD}] -nsteps $[4000*${BPD}]
#${FOLDER}/simulation_BPD${BPD}_FD -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/Herrmann_FD_${BPD} -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -bpdz ${BPD} -tend 1.25 -sim rti -fdump $[20*${BPD}] -nsteps $[4000*${BPD}]
done
#for BPD in 8 16 32
#do
#cp ../makefiles/simulation ${FOLDER}/simulation_BPD${BPD}_FD
#    export OMP_NUM_THREADS=48
#    bsub -n 48 -W 10:00 -o ${BASENAME}_BPD${BPD}_FD -J ${BASENAME}_BPD${BPD}_FD ${FOLDER}/simulation_BPD${BPD}_FD -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/Herrmann_FD_${BPD} -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -tend .9 -sim rti -fdump $[20*${BPD}] -nsteps $[2000*${BPD}]
#done


#cd ../makefiles/
#make clean
#make poisson=split-fftw bc=mixed config=production multiphase=true particles=true rk2=false -j
#cd -
#for BPD in 2 4
#do
#    cp ../makefiles/simulation ${FOLDER}/simulation_BPD${BPD}_P
#    export OMP_NUM_THREADS=48
#    bsub -n 48 -W 10:00 -o ${BASENAME}_BPD${BPD}_P -J ${BASENAME}_BPD${BPD}_P mpirun -np 1 ${FOLDER}/simulation_BPD${BPD}_P -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/Herrmann_P_${BPD} -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -tend .9 -sim rti -fdump $[20*${BPD}] -nsteps $[2000*${BPD}]
#done
#for BPD in 8 16
#do
#    cp ../makefiles/simulation ${FOLDER}/simulation_BPD${BPD}_P
#    export OMP_NUM_THREADS=48
#    bsub -n 48 -W 168:00 -o ${BASENAME}_BPD${BPD}_P -J ${BASENAME}_BPD${BPD}_P mpirun -np 32 ${FOLDER}/simulation_BPD${BPD}_P -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/Herrmann_P_${BPD} -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx ${BPD} -bpdy $[4*${BPD}] -tend .9 -sim rti -fdump $[20*${BPD}] -nsteps $[2000*${BPD}]
#done