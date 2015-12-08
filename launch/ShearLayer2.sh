cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=false densitydiff=true rk2=false -j;cd -
BASENAME2=ThinShearLayer3D_261115
FOLDER2=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}
mkdir ${FOLDER2}

for BPD in 2 4 8
#16 32 64
do
	cp ../makefiles/test ${FOLDER2}/ShearLayer_BPD${BPD}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME2}_BPD${BPD} -J ${BASENAME2}_BPD${BPD} ${FOLDER2}/ShearLayer_BPD${BPD} -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME2}/ShearLayer -CFL 0.1
done