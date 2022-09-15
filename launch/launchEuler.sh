#!/bin/bash
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
SETTINGSNAME=$1
BASENAME=$2
BASEPATH=${SCRATCH}/CubismUP3D/

source $SETTINGSNAME
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}
cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/simulation ${FOLDER}/simulation
cp -r ../source ${FOLDER}/

cd $FOLDER
unset LSB_AFFINITY_HOSTFILE
export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=1
echo $OPTIONS > settings.txt

if [ "${RUNLOCAL}" == "true" ] ; then
mpirun -n 360 --map-by core:PE=1 ./simulation ${OPTIONS} -factory-content "${FACTORY}"
else
bsub -J ${BASENAME} -W 24:00 -R "select[model==XeonGold_6150]fullnode" -n 360 "unset LSB_AFFINITY_HOSTFILE; mpirun -n 360 --map-by core:PE=1 ./simulation ${OPTIONS} -factory-content \"${FACTORY}\""
fi
