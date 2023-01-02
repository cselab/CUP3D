#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
SETTINGSNAME=$1
BASENAME=$2

WCLOCK=${WCLOCK:-24:00:00}
PARTITION=${PARTITION:-normal}
EXEC=${EXEC:-simulation}
BASEPATH="${SCRATCH}/CubismUP3D/"

source $SETTINGSNAME

if [ ${PSOLVER:0:4} == 'cuda' ] ; then
  export OMP_PLACES=cores
  export OMP_PROC_BIND=close
  export TASKS_PER_NODE=1
  export HINT=multithread
  if [ "${TASKS_PER_NODE}" -gt "1" ] ; then
    export CRAY_CUDA_MPS=1
  fi
  export OMP_NUM_THREADS=$(expr 12 / $TASKS_PER_NODE)
else
  export TASKS_PER_NODE=12
  export OMP_NUM_THREADS=1
  export HINT=nomultithread
fi

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/${EXEC} ${FOLDER}/simulation
cp -r ../source ${FOLDER}/
cp -r ../Cubism/include/Cubism ${FOLDER}/
cp $0 ${FOLDER}
git diff HEAD > ${FOLDER}/gitdiff.diff

cd ${FOLDER}

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=s1160
#SBATCH --job-name="${BASENAME}"
#SBATCH --time=${WCLOCK}
#SBATCH --partition=${PARTITION}
#SBATCH --constraint=gpu
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=${TASKS_PER_NODE}
#SBATCH --cpus-per-task=${OMP_NUM_THREADS}
#SBATCH --threads-per-core=1
#SBATCH --hint=${HINT}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

srun ./simulation ${OPTIONS} -factory-content $(printf "%q" "${FACTORY}")

EOF
chmod 755 daint_sbatch
sbatch daint_sbatch
