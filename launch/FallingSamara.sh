#module load gcc
export OMP_NUM_THREADS=48
export LD_LIBRARY_PATH=/cluster/home/infk/cconti/Visualization/Software/mitsuba-panos/mitsuba/dist:/cluster/apps/python/3.2.2/x86_64/lib64:/usr/local/brutus/lib64:/cluster/apps/python/3.2.2/x86_64/multithreadblas:/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib:/cluster/apps/gcc/4.8.2/lib64:/cluster/apps/fftw/3.3.3/x86_64/gcc_4.4.6/serial/lib64:/cluster/apps/boost/1.51.0/x86_64/serial/gcc_4.4.6/lib64:/cluster/apps/xerces/3.1.1/lib/:/cluster/home/infk/cconti/dependencies/glew-1.12.0/lib:/cluster/home/infk/cconti/usr/lib:/cluster/work/infk/hbabak/apps/hdf5-1.8.8_gcc_serial/lib/:/cluster/work/infk/cconti/VTK5.8_gcc/lib/vtk-5.8/:/cluster/work/infk/hbabak/numactl-2.0.7/:/cluster/work/infk/hbabak/apps/fftw-2.1.5/lib:/cluster/apps/blcr/0.8.4/x86_64/lib64:/cluster/apps/openmpi/1.6.2/x86_64/gcc_4.7.0/lib:/cluster/apps/lsf/8.0/linux2.6-glibc2.3-x86_64/lib:/cluster/work/infk/wvanrees/apps/TBB/tbb41_20120718oss/build/linux_intel64_gcc_cc4.7.0_libc2.12_kernel2.6.32_release/

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/
BASENAME=FallingSamara_0812_

cd ../makefiles/;make clean;make config=production bc=mixed precision=single particles=false dlm=true -j;cd ../launch/

for CFL in 0.1
do
for LAMBDA in 1
do
for B in 2 4 8 16 32
do
	OPTIONS=" -bpdx ${B} -bpdy $[4*${B}] -bpdz ${B} -CFL ${CFL} -uinf .01 -shape samara -rhoS 300. -fdump 250 -lambda ${LAMBDA} -dlm ${LAMBDA}"

	# Re 100
	NAME_RE100=Re100_DLM${LAMBDA}_CFL${CFL}_$[${B}*32]

	FOLDER=${BASEPATH}${BASENAME}${NAME_RE100}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	cp FallingSamara.sh ${FOLDER}
	bsub -n 48 -W 10:00 -o ${BASENAME}${NAME_RE100} -J ${BASENAME}${NAME_RE100} ${FOLDER}/simulation -file ${FOLDER}/${BASENAME}${NAME_RE100} -sim falling -Re 1100 -tend 1000 ${OPTIONS}
#	${FOLDER}/simulation -file ${FOLDER}/${BASENAME}${NAME_RE100} -sim falling -Re 1100 -tend 1000 ${OPTIONS}
done
done
done