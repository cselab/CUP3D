#module load gcc
P=$1

export OMP_NUM_THREADS=48

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/
BASENAME=FallingSamara_MPI_${P}11_2112_.001Inertia_

cd ../makefiles/;make clean;make CC=mpic++ config=production bc=mixed precision=double dlm=true nthreads=48 -j;cd ../launch/

for B in 2
#4
#4 8 16
do
for CFL in 0.01
#0.1
do
for LAMBDA in 1
do
for ISO in 0.0
#0.004 0.008 0.012
do
	OPTIONS=" -nprocsx ${P} -nprocsy 1 -nprocsz 1 -bpdx ${B} -bpdy $[4*${P}*${B}] -bpdz $[${P}*${B}] -CFL ${CFL} -shape samara -rhoS 250. -fdump 250 -lambda ${LAMBDA} -dlm ${LAMBDA} -nu 0.00001 -isosurface ${ISO}"

	NAME=DLM${LAMBDA}_CFL${CFL}_Iso${ISO}_$[${B}*32]

	FOLDER=${BASEPATH}${BASENAME}${NAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp FallingSamara.sh ${FOLDER}
cd ${FOLDER}
bsub -n $[48*${P}] -R 'span[ptile=48]' -W 10:00 -o ${BASENAME}${NAME} -J ${BASENAME}${NAME} mpirun -np ${P} -pernode ./simulation -file ${BASENAME}${NAME} -sim falling -tend 1000 ${OPTIONS}
#mpirun -np ${P} ./simulation -file ${BASENAME}${NAME} -sim falling -tend 1000 ${OPTIONS}
cd -
done
done
done
done