#!/bin/bash
##SBATCH --dependency=afterany:44961681
#SBATCH --job-name="render"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --account=s1160
#SBATCH --time=00:15:00
#SBATCH --constraint=gpu
#SBATCH --array=0-330:10
module load daint-gpu ParaView
srun -n $SLURM_NTASKS  --cpu_bind=sockets pvbatch script-movie.py --frames ${SLURM_ARRAY_TASK_ID}
