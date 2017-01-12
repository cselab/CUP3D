#!/bin/bash
#export OMP_SCHEDULE=dynamic
#export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple
export MYROUNDS=10000
export USEMAXTHREADS=1
#needed for thread safety:
#export MV2_ENABLE_AFFINITY=0
export LD_LIBRARY_PATH=/users/novatig/accfft/build_cuda/:$LD_LIBRARY_PATH
SETTINGSNAME=settings.sh
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

echo $OPTIONS > settings.txt

srun -n ${NNODE} --cpu_bind=none --ntasks-per-node=1 --cpus-per-task=12 ./simulation ${OPTIONS}
