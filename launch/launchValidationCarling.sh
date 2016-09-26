#!/bin/bash
#module load gcc
NNODE=$1

INTERACTIVE=0
if [ $# -gt 1 ] ; then
    if [ "${2}" = "node" ]; then
        echo "Running on current node"
        INTERACTIVE=1
    fi
fi

BASEPATH=/cluster/scratch_xp/public/novatig/CubismUP3D/
BASENAME=FlowPastCarlingFishRe0550_dump2d_0_
NPROCESSORS=$((${NNODE}*48))

CFL=0.5
LAMBDA=1e4
BPDX=20
BPDY=32
BPDZ=24

NAME_RUN=BPD${BPDX}_CFL${CFL}_${1}RANKS

FOLDER=${BASEPATH}${BASENAME}${NAME_RUN}
mkdir -p ${FOLDER}

cp factoryCarling ${FOLDER}/factory
cp ../makefiles/simulation ${FOLDER}
cp launchValidationCarling.sh ${FOLDER}
cp runValidationCarling.sh ${FOLDER}/run.sh

#CURRDIR=`pwd`
cd $FOLDER

if [ $INTERACTIVE -eq 1 ] ; then 
    export OMP_NUM_THREADS=48
    export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
    export PATH=/cluster/home/mavt/chatzidp/usr/mpich3/bin/:$PATH
    export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/fftw-3.3.5/lib/:$LD_LIBRARY_PATH
    #export LD_LIBRARY_PATH=/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.9.2/lib/:$LD_LIBRARY_PATH
    #export LD_LIBRARY_PATH=/cluster/apps/hdf5/1.8.12/x86_64/gcc_4.9.2/openmpi_1.6.5/lib/:$LD_LIBRARY_PATH
    #export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/fftw-3.3.5/build_with_openmpi/lib/:$LD_LIBRARY_PATH
    OPTIONS=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ} -nprocsx ${NNODE} -CFL ${CFL} -length 0.2 -lambda ${LAMBDA} -nu 0.00001"
    sort $LSB_DJOB_HOSTFILE | uniq  > lsf_hostfile      
    #mpirun -np ${NNODE} -bynode ./simulation -tend 8 ${OPTIONS} 
    mpich_run -np ${NNODE} -ppn 1 -bind-to none -launcher ssh -f lsf_hostfile valgrind ./simulation --tool=memcheck --track-origins=yes --leak-check=yes -tend 8 ${OPTIONS}
else
    bsub -n ${NPROCESSORS} -R span[ptile=48] -W 12:00 -J ${BASENAME}${NAME_RUN} ./run.sh $NNODE 
fi

   #cd $CURRDIR
   #valgrind --tool=memcheck --track-origins=yes --leak-check=yes
