for R in 8
do
for P in single
do
for SY in 1. 1.02
do
for SZ in 1.02
#0.98 1.
do
	cd ../makefiles/
	make clean
	make config=production nthreads=8 bc=mixed precision=${P} particles=false dlm=true hdf=true movingframe=true bs=16 bsx=${R} sy=${SY} sz=${SZ} -j
	mv simulation simulation_daint_${R}_${SY}_${SZ}_${P}
	cd ../launch/

	for CFL in 0.1
	do
	for DLM in 1
#.1 1
	do
	for RHO in 250
	do
		BPDYZ=16
		JOB=Samara_ParamStudy_512_${SY}_${SZ}_P${P}_CFL${CFL}_DLM${DLM}_rhoS${RHO}
		export SCR=/scratch/daint/cconti/${JOB}
		mkdir -p $SCR
		echo "#!/usr/bin/ksh
## submit with sbatch
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cconti@mavt.ethz.ch

#SBATCH --job-name=${JOB}
#SBATCH --account=s436
#SBATCH --partition=viz

cp ../makefiles/simulation_daint_${R}_${SY}_${SZ}_${P} $SCR/.
cd $SCR

source /opt/modules/default/init/bashexport

MPICH_MAX_THREAD_SAFETY=multiple
export MPICH_NEMESIS_ASYNC_PROGRESS=1
export OMP_WAIT_POLICY=PASSIVE

export OMP_NUM_THREADS=8
aprun -r 1 -n 32 -N 1 -d 8 ./simulation_daint_${R}_${SY}_${SZ}_${P} -file ${JOB} -sim falling -tend 1000 -nprocsx 32 -nprocsy 1 -nprocsz 1 -bpdx 1 -bpdy ${BPDYZ} -bpdz ${BPDYZ} -CFL ${CFL} -shape samara -rhoS ${RHO} -dlm ${DLM} -nu 0.000654041451953662 -isosurface .0 -tdump .5
" > submit_$((64*${R}))_${SY}_${SZ}_P${P}_CFL${CFL}_DLM${DLM}_rhoS${RHO}
		sbatch submit_$((64*${R}))_${SY}_${SZ}_P${P}_CFL${CFL}_DLM${DLM}_rhoS${RHO}
	done
	done
	done
done
done
done
done
