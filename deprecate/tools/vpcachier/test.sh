NAME=FallingSamaraFixed_long_DP_512_180216_Dora
export SCR=/scratch/daint/cconti/${NAME}

cd $SCR

for F in $(ls *.channel0)
do
	S=$(awk 'BEGIN {FS="[.]"} {print $1}' <<< $F)
	if [ ! -f $S.vpcache ]
	then
		if [[ $S=chi* ]]
		then
			echo $S
			#aprun -r 1 -n ${SLURM_NTASKS} -N ${SLURM_NTASKS_PER_NODE} -d 8 ./vpcachier -wtype_read 1 -wtype_write 3 -eps 0 -min 0.4 -max 0.6 -simdata $S -vp $S.vpcache
		fi
		if [[ $S=vorticity* ]]
		then
			echo $S
			#aprun -r 1 -n ${SLURM_NTASKS} -N ${SLURM_NTASKS_PER_NODE} -d 8 ./vpcachier -wtype_read 1 -wtype_write 3 -eps 0 -min 0 -max 500 -simdata $S -vp $S.vpcache
		fi
	fi
done
