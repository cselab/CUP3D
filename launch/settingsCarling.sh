#!/bin/bash
NNODEX=${NNODEX:-16}
NNODEY=${NNODEY:-1}
NNODE=$(($NNODEX * $NNODEY))

BPDX=${BPDX:-16}
BPDY=${BPDY:-2}
BPDZ=${BPDZ:-2}

NU=${NU:-0.00001}
BC=${BC:-freespace}

FACTORY='CarlingFish L=0.2 T=1.0 xpos=0.3 bFixToPlanar=1 bFixFrameOfRef=1 heightProfile=danio widthProfile=danio
'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 1.0 -tend 10.0 "
OPTIONS+=" -nslices 2 -slice1_direction 1 -slice2_direction 2 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.3 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 4 -levelStart 3 -Rtol 1.0 -Ctol 0.1"
OPTIONS+=" -Advection3rdOrder=true" 
