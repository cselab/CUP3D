BASENAME=ThinShearLayer_UP_180915_CFL_Euler_VarRho
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
mkdir ${FOLDER}

cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=false densitydiff=true rk2=false -j;cd -

for CFL in 0.05 0.1 0.3 0.4
do
    cp ../makefiles/test ${FOLDER}/test_FD
    export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_FD_CFL${CFL} -J ${BASENAME}_FD_CFL${CFL} ${FOLDER}/test_FD -minBPD 8 -maxBPD 8 -CFL ${CFL} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/ShearLayer_FD_CFL${CFL}
done

cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=true densitydiff=true rk2=false -j;cd -

for CFL in 0.05 0.1 0.3 0.4
do
cp ../makefiles/test ${FOLDER}/test_P
export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_P_CFL${CFL} -J ${BASENAME}_P_CFL${CFL} ${FOLDER}/test_P -minBPD 8 -maxBPD 8 -CFL ${CFL} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/ShearLayer_P_CFL${CFL}
done