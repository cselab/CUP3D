#!/bin/bash
if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
SETTINGSNAME=$1
BASENAME=$2

EXEC=${EXEC:-simulation}
BASEPATH="${SCRATCH}/CubismUP3D/"

source $SETTINGSNAME

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp ../bin/${EXEC} ${FOLDER}/simulation

cd ${FOLDER}

srun --ntasks-per-node=12 --nodes=$SLURM_NNODES --cpus-per-task=1 simulation ${OPTIONS} -factory-content $"${FACTORY}" | tee out.log
