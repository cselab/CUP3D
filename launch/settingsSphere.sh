#!/bin/bash
NNODE=64
DLM=${DLM:--1}
LAMBDA=${LAMBDA:-1e6}
BPDX=${BPDX:-16}
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
OPTIONS+=" -CFL 0.4 -lambda ${LAMBDA} -use-dlm ${DLM} -nu ${NU}"
OPTIONS+=" -ImplicitPenalization 1"
OPTIONS+=" -levelMax 4 -levelStart 2 -Rtol 0.40 -Ctol 0.10"
OPTIONS+=" -extentx 2.0 "
OPTIONS+=" -TimeOrder 2 "
