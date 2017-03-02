#!/bin/bash

MYNAME=`whoami`
BASENAME=test_static
NNODE=4
NPROCESSORS=$((${NNODE}*12))
FOLDER="./"
mkdir -p ./0
mkdir -p ./1

cat <<EOF >daint_sbatch
#!/bin/bash -l

#SBATCH --account=s658
#SBATCH --job-name="${BASENAME}"
#SBATCH --output=${BASENAME}_out_%j.txt
#SBATCH --error=${BASENAME}_err_%j.txt
#SBATCH --time=00:30:00
#SBATCH --nodes=${NNODE}
#SBATCH --partition=debug
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --threads-per-core=2
#SBATCH --constraint=gpu
export CRAY_CUDA_MPS=1


export LD_LIBRARY_PATH=/users/novatig/accfft/build_cuda/:$LD_LIBRARY_PATH
module load daint-gpu GSL/2.1-CrayGNU-2016.11 cray-hdf5-parallel/1.10.0
module load cudatoolkit/8.0.44_GA_2.2.7_g4a6c213-2.1 fftw/3.3.4.10
export OMP_NUM_THREADS=24
export MYROUNDS=10000
export USEMAXTHREADS=1
srun --ntasks ${NNODE} --threads-per-core=2 --ntasks-per-node=1 --cpus-per-task=12 ../ms
#./simulation ${OPTIONS}
EOF
#mpirun -l -n 8
chmod 755 daint_sbatch
sbatch daint_sbatch
