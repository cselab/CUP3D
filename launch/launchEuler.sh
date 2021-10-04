#!/bin/bash
SETTINGSNAME=$1
BASENAME=$2
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
BASEPATH=${SCRATCH}/CubismUP3D/

NTHREADS=1

if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

NPROCESSORS=$((${NNODE}*${NTHREADS}))
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/simulation ${FOLDER}/simulation

cd $FOLDER

unset LSB_AFFINITY_HOSTFILE
export OMP_NUM_THREADS=1 #${NTHREADS}
echo $OPTIONS > settings.txt

#bsub -J ${BASENAME} -W 24:00 -n 128 "unset LSB_AFFINITY_HOSTFILE; mpirun -n 128 --map-by core:PE=1 ./simulation ${OPTIONS} -factory-content \"${FACTORY}\""
bsub -J ${BASENAME} -W 24:00 -R "select[model==XeonGold_6150]fullnode" -n 360 "unset LSB_AFFINITY_HOSTFILE; mpirun -n 360 --map-by core:PE=1 ./simulation ${OPTIONS} -factory-content \"${FACTORY}\""

