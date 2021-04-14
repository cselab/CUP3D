#!/bin/bash
NNODE=64
DLM=${DLM:--1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-32} #16
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
LEVELS=${LEVELS:-4}
# Re=300 <-> NU=0.00005208333333; Re=500 <-> NU=0.00003125; Re=1000 <-> NU=0.000015625; Re=10'000 <-> NU=0.0000015625; Re=420'000 <-> NU=0.000000037202381; Re=1'140'000 <-> NU=0.00000001370614
NU=${NU:-0.00000001370614}
BC=${BC:-dirichlet}
FACTORY='Sphere L=0.125 xpos=0.6 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1'

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 1.0 -tend 100 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.25 -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax 4 -levelStart 3 -Rtol 0.40 -Ctol 0.10"
OPTIONS+=" -extentx 4.0 "
OPTIONS+=" -TimeOrder 2 "
OPTIONS+=" -poissonTol 1e-6 "

####
#### For a smaller domain use BPDX=16 and -extentx 2.0
####
