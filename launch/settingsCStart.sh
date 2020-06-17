#!/bin/bash
NNODEX=${NNODEX:-64}
NNODEY=${NNODEY:-1}
NNODE=$(($NNODEX * $NNODEY))

BPDX=${BPDX:-64}
BPDY=${BPDY:-$((${BPDX}))}
BPDZ=${BPDZ:-$((${BPDX/4}))}

NU=${NU:-0.0001931818182}
BC=${BC:-freespace}

FACTORY='CStartFish L=0.25 T=0.5882352941 xpos=0.50 bFixToPlanar=1 bFixFrameOfRef=0 heightProfile=cstart widthProfile=cstart
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.01 -tend 1.5882352941 "
#OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -use-dlm 10 -nu ${NU}"
