#!/bin/bash
#SBATCH --job-name="PV-samara"
#SBATCH --nodes=32
#SBATCH --ntasks=256
#SBATCH --partition=viz
#SBATCH --account=s436
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cconti@mavt.ethz.ch
#SBATCH --output=pvbatch-daint-OUT.log
#SBATCH --error=pvbatch-daint-ERR.log

#SBATCH --constraint=startx

export DISPLAY=:0
export LD_LIBRARY_PATH=/opt/cray/nvidia/default/lib64/:$LD_LIBRARY_PATH

for I in {0..5000}
do
	aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
done
