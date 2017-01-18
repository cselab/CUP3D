#!/bin/bash
srcpath=$1
destpath=$2

MYNAME=`whoami`
BASEPATH="/scratch/snx3000/${MYNAME}/CubismUP3D/"
#lfs setstripe -c 1 ${BASEPATH}${RUNFOLDER}
SETTINGSNAME=${srcpath}"/settings.sh"
if [ ! -f $SETTINGSNAME ];then
    echo ${SETTINGSNAME}" not found! - exiting"
    exit -1
fi
source $SETTINGSNAME
OPTIONS+=" -saveSplicer 1"
OPTIONS+=" -deserialization ${srcpath}"
OPTIONS+=" -serialization ${destpath}"

#ls ${srcpath}/restart_*status > list.dat
echo `cd ${srcpath} && ls restart_*status` > list.dat

OPTIONS+=" -list list.dat"
echo "running srun with options " ${OPTIONS}
sbatchf="daint_sbatch"

cat <<EOF >$sbatchf
#!/bin/bash -l

#SBATCH --account=s658
#SBATCH --time=12:00:00
#SBATCH --nodes=${NNODE}
# #SBATCH --partition=viz
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --constraint=gpu
export CRAY_CUDA_MPS=1
export LD_LIBRARY_PATH=/users/novatig/accfft/build_cuda/:$LD_LIBRARY_PATH
module load daint-gpu GSL/2.1-CrayGNU-2016.11 cray-hdf5-parallel/1.10.0
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1 fftw/3.3.4.10
export OMP_NUM_THREADS=24
export MYROUNDS=10000
export USEMAXTHREADS=1
srun -n ${NNODE} --ntasks-per-node=1 --cpus-per-task=24 /users/novatig/CubismUP_3D/makefiles/simulation ${OPTIONS}
EOF

chmod 755 $sbatchf
sbatch $sbatchf
