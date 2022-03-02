#!/bin/bash
NNODE=16

BPDX=${BPDX:-8}
BPDY=${BPDY:-4}
BPDZ=${BPDZ:-4}

NU=${NU:-0.00001}
BC=${BC:-freespace}

FACTORY='CarlingFish L=0.2 T=1.0 xpos=0.6 planarAngle=0.0 bFixFrameOfRef=1 heightProfile=stefan widthProfile=stefan
'
FACTORY+='CarlingFish L=0.2 T=1.0 xpos=0.6 ypos=0.15 planarAngle=-25.0 heightProfile=stefan widthProfile=stefan
'

OPTIONS=
OPTIONS+=" -extentx 1.0"
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -tdump 0.05 -tend 20.0"
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.4 -nu ${NU}"
OPTIONS+=" -levelMax 4 -levelStart 2 -Rtol 0.5 -Ctol 0.05"
