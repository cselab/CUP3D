#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Usage "$0" SETTINGSNAME BASENAME"
  exit 1
fi
SETTINGSNAME=$1
BASENAME=$2

WCLOCK=${WCLOCK:-48:00:00}
PARTITION=${PARTITION:-standard}
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
  export OMP_NUM_THREADS=$(expr 128 / $TASKS_PER_NODE)
else
  export TASKS_PER_NODE=128
  export OMP_NUM_THREADS=1
fi

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/${EXEC} ${FOLDER}/simulation
cp -r ../source ${FOLDER}/
cp -r ../Cubism/include/Cubism ${FOLDER}/
git diff HEAD > ${FOLDER}/gitdiff.diff

cd ${FOLDER}

cat <<EOF >lumi_sbatch
#!/bin/bash -l

#SBATCH --account=${ACCOUNT}
#SBATCH --job-name="${BASENAME}"
#SBATCH --time=${WCLOCK}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=${TASKS_PER_NODE}
#SBATCH --cpus-per-task=${OMP_NUM_THREADS}
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
export FI_CXI_RX_MATCH_MODE=hybrid
srun ./simulation ${OPTIONS} -factory-content $(printf "%q" "${FACTORY}")

EOF
chmod 755 lumi_sbatch
sbatch lumi_sbatch
