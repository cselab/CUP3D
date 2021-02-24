#!/bin/bash
NNODE=128

BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}

NU=${NU:-0.00000004} # Re = 1,000,000
BC=${BC:-freespace}

FACTORY='StefanFish L=0.2 T=1.0 xpos=0.3 bFixFrameOfRef=1 heightProfile=danio widthProfile=danio Correct=true
'

OPTIONS=
OPTIONS+=" -extentx 1.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.25 -tend 50.0 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.5 -use-dlm 10 -nu ${NU}"
OPTIONS+=" -levelMax 5 -levelStart 3 -Rtol 0.5 -Ctol 0.05"
OPTIONS+=" -implicitPenalization 1"
OPTIONS+=" -TimeOrder 2"
