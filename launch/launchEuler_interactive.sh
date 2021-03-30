#!/bin/bash
SETTINGSNAME=$1
BASENAME=$2
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
BASEPATH=${SCRATCH}/RUNS/
#BASEPATH=~/RUNS/

INTERACTIVE=0
NNODE=1

if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/simulation ${FOLDER}/simulation

cd $FOLDER
unset LSB_AFFINITY_HOSTFILE
export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=36
echo $OPTIONS > settings.txt
export LD_LIBRARY_PATH=/cluster/home/novatig/hdf5-1.10.1/gcc_6.3.0_openmpi_2.1/lib/:$LD_LIBRARY_PATH

mpirun -n 5 --map-by node:PE=36 ./simulation ${OPTIONS} -factory-content "${FACTORY}"
