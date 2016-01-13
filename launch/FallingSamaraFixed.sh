#module load gcc
P=$1

export OMP_NUM_THREADS=48
#export LD_LIBRARY_PATH=/cluster/home/infk/cconti/Visualization/Software/mitsuba-panos/mitsuba/dist:/cluster/apps/python/3.2.2/x86_64/lib64:/usr/local/brutus/lib64:/cluster/apps/python/3.2.2/x86_64/multithreadblas:/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib:/cluster/apps/gcc/4.8.2/lib64:/cluster/apps/fftw/3.3.3/x86_64/gcc_4.4.6/serial/lib64:/cluster/apps/boost/1.51.0/x86_64/serial/gcc_4.4.6/lib64:/cluster/apps/xerces/3.1.1/lib/:/cluster/home/infk/cconti/dependencies/glew-1.12.0/lib:/cluster/home/infk/cconti/usr/lib:/cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc_serial/lib/:/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:/cluster/work/infk/hbabak/numactl-2.0.7/:/cluster/work/infk/hbabak/apps/fftw-2.1.5/lib:/cluster/apps/blcr/0.8.4/x86_64/lib64:/cluster/apps/openmpi/1.6.2/x86_64/gcc_4.7.0/lib:/cluster/apps/lsf/8.0/linux2.6-glibc2.3-x86_64/lib:/cluster/work/infk/wvanrees/apps/TBB/tbb41_20120718oss/build/linux_intel64_gcc_cc4.7.0_libc2.12_kernel2.6.32_release/

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/
BASENAME=FallingSamaraFixed_tdump_MPI${P}_1201_1Inertia_Wang_testCoMPos

cd ../makefiles/
make clean
make CC=mpic++ config=production bc=mixed precision=double particles=false dlm=true hdf=true movingframe=true -j
cd ../launch/

for B in 1
#2 4 8
do
for CFL in 0.1
#0.01 0.1
do
for LAMBDA in 1.
#.1 1.
do
for ISO in 0.004
#0.0 0.004
do
for NU in 0.0000001
#0.001 0.00001 0.0000001
do
	OPTIONS=" -nprocsx ${P} -nprocsy 1 -nprocsz 1 -bpdx ${B} -bpdy $[${P}*${B}] -bpdz $[${P}*${B}] -CFL ${CFL} -shape samara -rhoS 250. -lambda ${LAMBDA} -dlm ${LAMBDA} -nu ${NU} -isosurface ${ISO} -tdump .01 -fdump 500"

	NAME=DLM${LAMBDA}_CFL${CFL}_Iso${ISO}_$[${P}*${B}*32]_Nu${NU}

	FOLDER=${BASEPATH}${BASENAME}${NAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	cp FallingSamaraFixed.sh ${FOLDER}
	cd ${FOLDER}
	bsub -n $[48*${P}] -R 'span[ptile=48]' -W 40:00 -o ${BASENAME}${NAME} -J ${BASENAME}${NAME} mpirun -np ${P} -pernode ./simulation -file ${BASENAME}${NAME} -sim falling -tend 1000 ${OPTIONS}
#	./simulation -file ${BASENAME}${NAME} -sim falling -tend 1000 ${OPTIONS}
	cd -
done
done
done
done
done