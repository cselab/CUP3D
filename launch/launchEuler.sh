#!/bin/bash
SETTINGSNAME=$1
BASEPATH=${SCRATCH}/CubismUP3D/

INTERACTIVE=0
#INTERACTIVE=1
if [ $# -gt 1 ] ; then
    if [ "${2}" = "node" ]; then
        echo "Running on current node"
        INTERACTIVE=1
    fi
fi

if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

NPROCESSORS=$((${NNODE}*12))
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
cp ${FFACTORY} ${FOLDER}/factory
cp ../makefiles/simulation ${FOLDER}
#cp launchBrutus.sh ${FOLDER}
cp runEuler.sh ${FOLDER}/run.sh

#CURRDIR=`pwd`
cd $FOLDER

unset LSB_AFFINITY_HOSTFILE
export OMP_NUM_THREADS=24
export MV2_ENABLE_AFFINITY=0
echo $OPTIONS > settings.txt

export LD_LIBRARY_PATH=/cluster/home/novatig/hdf5-1.10.1/gcc_6.3.0_openmpi_2.1/lib/:$LD_LIBRARY_PATH

if [ $INTERACTIVE -eq 1 ] ; then 
   #mpirun -n ${NNODE} --map-by ppr:1:socket:pe=12 --bind-to core -report-bindings --mca mpi_cuda_support 0 valgrind --tool=memcheck --leak-check=yes --track-origins=yes --show-reachable=yes ./simulation ${OPTIONS}
   #mpirun -n ${NNODE} --map-by ppr:1:socket:pe=12 --bind-to core -report-bindings --mca mpi_cuda_support 0 valgrind --tool=memcheck --undef-value-errors=no --num-callers=500  ./simulation ${OPTIONS}
   mpirun -n ${NNODE} --map-by ppr:1:socket:pe=12 --bind-to core -report-bindings --mca mpi_cuda_support 0  ./simulation ${OPTIONS}
   #mpirun -np ${NNODE} -ppn 1 ./simulation ${OPTIONS}
else
    bsub  -R "rusage[mem=320] select[model==XeonE5_2680v3]" -n ${NPROCESSORS} -W 24:00 -J ${BASENAME} < run.sh
fi

#cd $CURRDIR
#valgrind --tool=memcheck --track-origins=yes --leak-check=yes
