
function findvort
{
    tstep=$1
    P1="/scratch/daint/cconti/CVT/CVT_Re10k_M.1/cached/vorticity${tstep}.vpcache"

    if [ -f "$P1" ]
    then
	echo $P1
    fi
}


for P in $(ls /scratch/daint/cconti/CVT/CVT_Re10k_M.1/cached/*.vpcache)
do
    F=$(basename $P)
    F1=${F%.vpcache}
    F2=${F1#vorticity}
    echo $F2
done | sort -n -k1 | uniq > /tmp/asdo.txt 

for T in $(cat /tmp/asdo.txt)
do
    findvort $T
done
