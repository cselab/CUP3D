#!/bin/bash
#module load gcc
NNODE=$1
NTHREADSPERNODE=48

INTERACTIVE=0
if [ $# -gt 1 ] ; then
    if [ "${2}" = "node" ]; then
        echo "Running on current node"
        INTERACTIVE=1
    fi
fi

BASEPATH=/cluster/scratch_xp/public/novatig/CubismUP3D/
BASENAME=FlowPastFixedSphereRe200_Validate3_
NPROCESSORS=$((${NNODE}*${NTHREADSPERNODE}))

CFL=0.1
LAMBDA=1e4
BPDX=20
BPDY=32

NAME_RUN=BPD${BPDX}_CFL${CFL}

FOLDER=${BASEPATH}${BASENAME}${NAME_RUN}
mkdir -p ${FOLDER}

cp factoryFixedSphere ${FOLDER}/factory
cp ../makefiles/simulation ${FOLDER}
cp launchValidationSphere.sh ${FOLDER}
cp runValidationSphere.sh ${FOLDER}/run.sh
#CURRDIR=`pwd`
cd $FOLDER

if [ $INTERACTIVE -eq 1 ] ; then 
    export OMP_NUM_THREADS=$NTHREADSPERNODE
    export LD_LIBRARY_PATH=/cluster/home/mavt/chatzidp/usr/mpich3/lib/:$LD_LIBRARY_PATH
    export PATH=/cluster/home/mavt/chatzidp/usr/mpich3/bin/:$PATH
    export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/fftw-3.3.5/lib/:$LD_LIBRARY_PATH
    OPTIONS=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDY} -nprocsx ${NNODE} -CFL ${CFL} -length 0.1 -uinfx 0.1 -lambda ${LAMBDA} -nu 0.00005"
    sort $LSB_DJOB_HOSTFILE | uniq  > lsf_hostfile      
    mpich_run -np ${NNODE} -ppn 1 -bind-to none -launcher ssh -f lsf_hostfile ./simulation -tend 8 ${OPTIONS}
else
    bsub -n ${NPROCESSORS} -R span[ptile=48] -sp 100 -W 12:00 -J ${BASENAME}${NAME_RUN} ./run.sh $NNODE 
fi

#cd $CURRDIR

    #valgrind --tool=memcheck --track-origins=yes --leak-check=yes 

    
