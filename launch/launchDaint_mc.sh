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

FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/${EXEC} ${FOLDER}/simulation
cp -r ../source ${FOLDER}/
cp -r ../Cubism ${FOLDER}/
cp $0 ${FOLDER}
git diff HEAD > ${FOLDER}/gitdiff.diff

cd ${FOLDER}

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=eth2
#SBATCH --job-name="${BASENAME}"
#SBATCH --output=${BASENAME}_out_%j.txt
#SBATCH --error=${BASENAME}_err_%j.txt
#SBATCH --time=${WCLOCK}
#SBATCH --partition=${PARTITION}
#SBATCH --constraint=mc
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=18
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=18

srun ./simulation ${OPTIONS} -factory-content $(printf "%q" "${FACTORY}")

EOF
chmod 755 daint_sbatch
sbatch daint_sbatch
