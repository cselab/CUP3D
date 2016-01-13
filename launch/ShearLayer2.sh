P=$1

cd ../makefiles/;make clean;make CC=mpic++ bc=periodic config=production multiphase=true precision=double vertexcentered=true poisson=split-fftw particles=false densitydiff=true rk2=false -j;cd -
BASENAME2=ThinShearLayer3D_120116_MPI2
FOLDER2=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}
mkdir ${FOLDER2}

for BPD in 4 8 16
#32 64
do
cp ../makefiles/test ${FOLDER2}/ShearLayer_BPD${BPD}
cd ${FOLDER2}
export OMP_NUM_THREADS=48;bsub -n $[48*${P}] -R 'span[ptile=48]' -W 10:00 -o ${BASENAME2}_BPD${BPD} -J ${BASENAME2}_BPD${BPD} mpirun -np ${P} -pernode ./ShearLayer_BPD${BPD} -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file ShearLayer -CFL 0.1
#mpirun -np 1 -pernode ./ShearLayer_BPD${BPD} -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file ShearLayer -CFL 0.1
cd -
done