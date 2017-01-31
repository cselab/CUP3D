#!/bin/bash
SETTINGSNAME=$1
BASENAME=$2
if [ $# -gt 3 ] ; then
  NNODEX=$3
  NNODEY=$4
  NNODE=$((${3}*${4}))
  BASENAME=${BASENAME}_NX${3}_NY${4}
else
  if [ $# -gt 2 ] ; then
    NNODE=$3
    NNODEX=$NNODE
    NNODEY=1
  else
    echo "Specify number of nodes!"
    exit -1
  fi
fi
echo $BASENAME $NNODE

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
#cp ../makefiles/simulation ${FOLDER}/simulation
#cp ../makefiles/simulation+pat ${FOLDER}/simulation
cp simulation+apa ${FOLDER}/simulation
cp launchDaint.sh ${FOLDER}
cp runDaint.sh ${FOLDER}/run.sh

cd ${FOLDER}

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=s658
#SBATCH --job-name="${BASENAME}"
#SBATCH --output=${BASENAME}_%j.txt
#SBATCH --error=${BASENAME}_%j.txt
#SBATCH --time=01:00:00
#SBATCH --nodes=${NNODE}
# #SBATCH --partition=viz
#SBATCH --threads-per-core=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu
export CRAY_CUDA_MPS=1
# #module load perftools-base/6.4.3 perftools

export LD_LIBRARY_PATH=/users/novatig/accfft/build_shared/:$LD_LIBRARY_PATH
module load daint-gpu GSL/2.1-CrayGNU-2016.11 cray-hdf5-parallel/1.10.0
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1 fftw/3.3.4.10
export OMP_NUM_THREADS=24
export MYROUNDS=10000
export USEMAXTHREADS=1
srun --ntasks ${NNODE} --threads-per-core=2 --ntasks-per-node=1 --cpus-per-task=12 ./simulation ${OPTIONS}
EOF

chmod 755 daint_sbatch
sbatch daint_sbatch
