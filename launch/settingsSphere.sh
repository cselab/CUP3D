#!/bin/bash
NNODE=64
DLM=${DLM:--1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-32}
BPDY=${BPDY:-8}
BPDZ=${BPDZ:-8}
NU=${NU:-0.00005208333333} # Re 300, to compare against Mimeau Cottet Mortazavi 2016
BC=${BC:-dirichlet}
FACTORY='Sphere L=0.125 xpos=0.6 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1
'
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
OPTIONS+=" -poissonTolRel 1e-4 "

####
#### For a smaller domain use BPDX=16 and -extentx 2.0
####
