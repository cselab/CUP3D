#!/bin/bash
SETTINGSNAME=$1

WCLOCK=${WCLOCK:-24:00:00}
PARTITION=${PARTITION:-normal}
EXEC=${EXEC:-simulation}

MYNAME=`whoami`
BASEPATH="${SCRATCH}/CubismUP3D/"
#lfs setstripe -c 1 ${BASEPATH}${RUNFOLDER}

if [ ! -f $SETTINGSNAME ]; then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

NPROCESSORS=$((${NNODE}*12))
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
[[ -n "${FFACTORY}" ]] && cp ${FFACTORY} ${FOLDER}/factory
cp ../bin/${EXEC} ${FOLDER}/simulation
cp -r ../source ${FOLDER}/
cp $0 ${FOLDER}

git diff HEAD > ${FOLDER}/gitdiff.diff

cd ${FOLDER}

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=s658
#SBATCH --job-name="${BASENAME}"
#SBATCH --output=${BASENAME}_out_%j.txt
#SBATCH --error=${BASENAME}_err_%j.txt

#SBATCH --time=${WCLOCK}
#SBATCH --partition=${PARTITION}

#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --threads-per-core=1
#SBATCH --constraint=gpu
#SBATCH --mail-user=${MYNAME}@ethz.ch
#SBATCH --mail-type=ALL

module load daint-gpu GSL cray-hdf5-parallel cray-fftw
module load cudatoolkit/9.0.103_3.7-6.0.4.1_2.1__g72b395b

export MPICH_MAX_THREAD_SAFETY=multiple
export OMP_NUM_THREADS=12
srun --ntasks ${NNODE} --threads-per-core=1 --ntasks-per-node=1 --cpus-per-task=12 time ./simulation ${OPTIONS} -factory-content $(printf "%q" "${FACTORY}")
EOF

chmod 755 daint_sbatch
sbatch daint_sbatch
