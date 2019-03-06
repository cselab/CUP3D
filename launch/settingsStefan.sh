#!/bin/bash
NNODEX=${NNODEX:-4}
NNODEY=${NNODEY:-1}
NNODE=$(($NNODEX * $NNODEY))
DLM=${DLM:-1}
LAMBDA=${LAMBDA:-1e4}
BPDX=${BPDX:-64}
BPDY=${BPDY:-$((${BPDX}/2))} #${BPDY:-32}
BPDZ=${BPDZ:-$((${BPDX}/4))} #${BPDZ:-32}

NU=${NU:-0.0001636363636} # Re 4k
BC=${BC:-freespace}

#FACTORY='CarlingFish L=0.3 T=1 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1
#FACTORY='StefanFish L=0.3 T=1 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 Correct=1
FACTORY='CarlingFish L=0.25 T=1.0 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 1 -dump3D 1 -tdump 1 -tend 5 "
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
