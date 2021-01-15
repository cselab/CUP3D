#!/bin/bash
SETTINGSNAME=$1
BASENAME=$2
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
BASEPATH=${SCRATCH}/CubismUP3D/

NTHREADS=36

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
export OMP_NUM_THREADS=${NTHREADS}
echo $OPTIONS > settings.txt

bsub -J ${BASENAME} -R "select[model==XeonGold_6150]fullnode" -n ${NPROCESSORS} "unset LSB_AFFINITY_HOSTFILE; mpirun -n ${NNODE} --map-by node:PE=36 ./simulation ${OPTIONS} -factory-content \"${FACTORY}\""

