#!/bin/bash
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
SETTINGSNAME=$1
BASENAME=$2
BASEPATH=${SCRATCH}/RUNS/

source $SETTINGSNAME
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}
cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/simulation ${FOLDER}/simulation

cd $FOLDER
unset LSB_AFFINITY_HOSTFILE
export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=1
echo $OPTIONS > settings.txt
export LD_LIBRARY_PATH=/cluster/home/novatig/hdf5-1.10.1/gcc_6.3.0_openmpi_2.1/lib/:$LD_LIBRARY_PATH

mpirun -n 36 --map-by core:PE=1 ./simulation ${OPTIONS} -factory-content "${FACTORY}"
