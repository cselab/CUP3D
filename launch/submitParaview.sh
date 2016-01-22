#!/bin/bash
#SBATCH --job-name="paraview_samara"
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=8
#SBATCH --partition=viz
#SBATCH --account=s436
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cconti@mavt.ethz.ch
#SBATCH --output=/users/cconti/CubismUP_3D/launch/pvbatch-daint2-OUT.log
#SBATCH --error=/users/cconti/CubismUP_3D/launch/pvbatch-daint2-ERR.log

#SBATCH --constraint=startx

export DISPLAY=:0
export LD_LIBRARY_PATH=/opt/cray/nvidia/default/lib64/:$LD_LIBRARY_PATH

# both versions passed on nov 13, 2015, except with ghost-cell artifacts at processor boundaries
aprun -n $SLURM_NTASKS -N ${SLURM_NTASKS_PER_NODE} -d 1 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateMovies.py
