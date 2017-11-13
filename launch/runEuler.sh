#!/bin/bash
export OMP_NUM_THREADS=12
#export OMP_SCHEDULE=dynamic
#export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple
export MYROUNDS=10000
export USEMAXTHREADS=1
#needed for thread safety:
export MV2_ENABLE_AFFINITY=0
unset LSB_AFFINITY_HOSTFILE

SETTINGSNAME=settings.sh
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME
echo ${NNODE}
echo $OPTIONS > settings.txt

mpirun -n ${NNODE} --map-by ppr:1:socket:pe=12 --bind-to core -report-bindings --mca mpi_cuda_support 0 ./simulation ${OPTIONS}
