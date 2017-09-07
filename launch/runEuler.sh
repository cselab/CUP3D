#!/bin/bash
export OMP_NUM_THREADS=48
#export OMP_SCHEDULE=dynamic
#export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple
export MYROUNDS=10000
export USEMAXTHREADS=1
#needed for thread safety:
export MV2_ENABLE_AFFINITY=0

SETTINGSNAME=settings.sh
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME
echo ${NNODE}
echo $OPTIONS > settings.txt

mpirun -np ${NNODE} -ppn 1 ./simulation ${OPTIONS}
