#!/bin/bash
export OMP_NUM_THREADS=24
#needed for thread safety:
export MV2_ENABLE_AFFINITY=0

SETTINGSNAME=settings.sh
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

echo $OPTIONS > settings.txt

mpirun -np ${NNODE} -ppn 1 ./simulation ${OPTIONS}
