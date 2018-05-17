
function findchi
{
    tstep=$1
    P1="/scratch/daint/cconti/FallingSamaraFixed_long_DP_512_180216_Dora/cached/chi${tstep}.vpcache"
    P2="/scratch/daint/cconti/FallingSamaraFixed_long_DP_512_180216_Dora/cached/vorticity${tstep}.vpcache"

    if [ -f "$P1" ] && [ -f "$P2" ]
    then
	echo $P2
	echo $P1
    fi
}


for P in $(ls /scratch/daint/cconti/FallingSamaraFixed_long_DP_512_180216_Dora/cached/*.vpcache)
do
    F=$(basename $P)
    F1=${F%.vpcache}
    F2=${F1#vorticity}
    F3=${F2#chi}
    echo $F3
done | sort -n -k1 | uniq > /tmp/asdo.txt 

for T in $(cat /tmp/asdo.txt)
do
    findchi $T
done
