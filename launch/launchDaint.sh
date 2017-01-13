#!/bin/bash
SETTINGSNAME=$1

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
#SBATCH --output=${BASENAME}_%j.txt
#SBATCH --error=${BASENAME}_%j.txt
#SBATCH --time=00:30:00
#SBATCH --nodes=${NNODE}
# #SBATCH --partition=viz
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=gpu
export CRAY_CUDA_MPS=1
module swap PrgEnv-cray PrgEnv-gnu
module load daint-gpu fftw/3.3.4.10 cray-hdf5-parallel/1.10.0 cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1
export OMP_NUM_THREADS=24
export MYROUNDS=10000
export USEMAXTHREADS=1
#needed for thread safety:
#export MV2_ENABLE_AFFINITY=0
export LD_LIBRARY_PATH=/users/novatig/accfft/build_cuda/:$LD_LIBRARY_PATH
if [ ! -f settings.sh ];then
    echo "settings.sh not found! - exiting"
    exit -1
fi
source settings.sh
echo $OPTIONS > settings.txt
srun -n ${NNODE} --ntasks-per-node=1 --cpus-per-task=24 ./simulation ${OPTIONS}
EOF

chmod 755 daint_sbatch

sbatch daint_sbatch
