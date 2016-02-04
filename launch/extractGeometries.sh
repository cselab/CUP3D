#!/bin/bash
#SBATCH --job-name="PV-geo"
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --partition=viz
#SBATCH --account=s436
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cconti@mavt.ethz.ch
#SBATCH --output=pvbatch-daint-OUT.log
#SBATCH --error=pvbatch-daint-ERR.log

#SBATCH --constraint=startx

export DISPLAY=:0
export LD_LIBRARY_PATH=/opt/cray/nvidia/default/lib64/:$LD_LIBRARY_PATH

aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py

aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py

aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py
aprun -n $SLURM_NTASKS -N 8 `which pvbatch` --disable-xdisplay-test /users/cconti/CubismUP_3D/launch/generateGeometries.py

