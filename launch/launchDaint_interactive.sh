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

if [ ${PSOLVER:0:4} == 'cuda' ] ; then
  export OMP_PLACES=cores
  export OMP_PROC_BIND=close
  export TASKS_PER_NODE=1
  if [ "${TASKS_PER_NODE}" -gt "1" ] ; then
    export CRAY_CUDA_MPS=1
  fi
  export OMP_NUM_THREADS=$(expr 12 / $TASKS_PER_NODE)
else
  export TASKS_PER_NODE=12
  export OMP_NUM_THREADS=1
fi


FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp ../bin/${EXEC} ${FOLDER}/simulation

cd ${FOLDER}

srun --ntasks-per-node=$TASKS_PER_NODE --nodes=$SLURM_NNODES --cpus-per-task=$OMP_NUM_THREADS simulation ${OPTIONS} -factory-content $"${FACTORY}" | tee -a out.log
