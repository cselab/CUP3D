#!/bin/bash
NNODE=64
DLM=${DLM:-0}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-8}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
NU=${NU:-0.00005208333333} # Re 300, to compare against Mimeau Cottet Mortazavi 2016:
BC=${BC:-dirichlet}
FACTORY='Sphere L=0.125 xpos=0.125 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1
'
OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 1.0 -tend 12000 "
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -CFL 0.1 -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -Advection3rdOrder=true"
OPTIONS+=" -levelMax 4 -levelStart 3 -Rtol 0.40 -Ctol 0.10"
