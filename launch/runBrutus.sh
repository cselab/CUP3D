#!/bin/bash
export OMP_NUM_THREADS=48
#export OMP_SCHEDULE=dynamic
#export MPICH_NEMESIS_ASYNC_PROGRESS=1
#export MPICH_MAX_THREAD_SAFETY=multiple
export MYROUNDS=10000
export USEMAXTHREADS=1
#needed for thread safety:
export MV2_ENABLE_AFFINITY=0
#modules compiled for gcc 5.2.0
export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/hdf5-1.8.17/mvapich2.2.1/${_MY_CC_}/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/fftw-3.3.5/mvapich2.2.1/${_MY_CC_}/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cluster/home03/mavt/novatig/accfft/mvapich2.2.1/${_MY_CC_}/:$LD_LIBRARY_PATH
SETTINGSNAME=settings.sh
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

echo $OPTIONS > settings.txt

mpirun -np ${NNODE} -ppn 1  ./simulation ${OPTIONS}
