#cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=fftw particles=false densitydiff=false rk2=true -j;cd -
#BASENAME=ThinShearLayer_140915_RK2_ConstRho
#FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
#mkdir ${FOLDER}

#for BPD in 8 16 32
#do
#   cp ../makefiles/test ${FOLDER}/test_BPD${BPD}_FD
#    export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 -o ${BASENAME}_BPD${BPD}_FD -J ${BASENAME}_BPD${BPD}_FD ${FOLDER}/test_BPD${BPD}_FD -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/ShearLayer_CFL0.05_FD
#done

#cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=fftw particles=true densitydiff=false rk2=true -j;cd -

#for BPD in 8 16 32
#do
#   cp ../makefiles/test ${FOLDER}/test_BPD${BPD}_P
#   export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 -o ${BASENAME}_BPD${BPD}_FD -J ${BASENAME}_BPD${BPD}_FD ${FOLDER}/test_BPD${BPD}_FD -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/ShearLayer_CFL0.05_P
#done


cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=false densitydiff=true rk2=false -j;cd -
BASENAME2=ThinShearLayer_UP_180915_nuRho
FOLDER2=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}
mkdir ${FOLDER2}

for BPD in 8 16 32 64
do
	cp ../makefiles/test ${FOLDER2}/test_BPD${BPD}_FD_nu5
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME2}_BPD${BPD}_FD_nu5 -J ${BASENAME2}_BPD${BPD}_FD_nu5 ${FOLDER2}/test_BPD${BPD}_FD_nu5 -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}/ShearLayer_FD_nu5 -CFL 0.5
done

cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=true densitydiff=true rk2=false -j;cd -

for BPD in 8 16 32 64
do
	cp ../makefiles/test ${FOLDER2}/test_BPD${BPD}_particles_nu5
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME2}_BPD${BPD}_particles_nu5 -J ${BASENAME2}_BPD${BPD}_particles_nu5 ${FOLDER2}/test_BPD${BPD}_particles_nu5 -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}/ShearLayer_particles_nu5 -CFL 0.5
done