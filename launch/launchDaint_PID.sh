#!/bin/bash
SETTINGSNAME=$1
FACTORYNAME=$2

MYNAME=`whoami`
BASEPATH="/scratch/snx3000/${MYNAME}/CubismUP3D/"
#lfs setstripe -c 1 ${BASEPATH}${RUNFOLDER}

if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME

NPROCESSORS=$((${NNODE}*12))
FOLDER=${BASEPATH}${BASENAME}
mkdir -p ${FOLDER}

cp $SETTINGSNAME ${FOLDER}/settings.sh
cp ${FFACTORY} ${FOLDER}/factory
cp ../makefiles/simulation ${FOLDER}
cp launchDaint.sh ${FOLDER}
cp runDaint.sh ${FOLDER}/run.sh

cd ${FOLDER}

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=s658
#SBATCH --job-name="${BASENAME}"
#SBATCH --output=${BASENAME}_out_%j.txt
#SBATCH --error=${BASENAME}_err_%j.txt
#SBATCH --time=24:00:00
#SBATCH --nodes=${NNODE}
# #SBATCH --partition=debug
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --threads-per-core=2
#SBATCH --constraint=gpu
export CRAY_CUDA_MPS=1


export LD_LIBRARY_PATH=/users/novatig/accfft/build_dbg/:$LD_LIBRARY_PATH
module load daint-gpu GSL/2.1-CrayGNU-2016.11 cray-hdf5-parallel/1.10.0
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1 fftw/3.3.4.10
export OMP_NUM_THREADS=24
export MYROUNDS=10000
export USEMAXTHREADS=1
srun --ntasks ${NNODE} --threads-per-core=2 --ntasks-per-node=1 --cpus-per-task=12 ./simulation ${OPTIONS}
#srun --ntasks ${NNODE} --threads-per-core=2 --ntasks-per-node=1 --cpus-per-task=12 valgrind --trace-syscalls=yes --num-callers=100 --tool=memcheck --leak-check=yes --track-origins=yes --show-reachable=yes ./simulation ${OPTIONS}
EOF

chmod 755 daint_sbatch
sbatch daint_sbatch
