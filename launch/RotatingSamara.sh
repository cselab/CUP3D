#module load gcc
export OMP_NUM_THREADS=48

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/
BASENAME=RotatingSamara_140116_

cd ../makefiles/;make clean;make CC=mpic++ config=production bc=periodic precision=double particles=false dlm=true hdf=true -j;cd ../launch/

for B in 8
#4 8 16
do
for CFL in 0.01
#0.1
do
for LAMBDA in 1
do
for ISO in 0.01
#0.004 0.008 0.012
do
	OPTIONS="-nprocsx 1 -nprocsy 1 -nprocsz 1 -bpdx ${B} -bpdy ${B} -bpdz ${B} -CFL ${CFL} -shape samara -fdump 100 -lambda ${LAMBDA} -dlm ${LAMBDA} -nu 0.000042 -isosurface ${ISO}"

	NAME=DLM${LAMBDA}_CFL${CFL}_Iso${ISO}_$[${B}*32]

	FOLDER=${BASEPATH}${BASENAME}${NAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	cp RotatingSamara.sh ${FOLDER}
	cd ${FOLDER}
#	bsub -n 48 -W 10:00 -o ${BASENAME}${NAME} -J ${BASENAME}${NAME} ./simulation -file ${BASENAME}${NAME} -sim falling -tend 1000 ${OPTIONS}
	mpirun -np 1 ./simulation -file ${BASENAME}${NAME} -sim rotating -tend 1000 ${OPTIONS}
	cd -
done
done
done
done