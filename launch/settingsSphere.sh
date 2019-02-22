#!/bin/bash
NNODEX=${NNODEX:-8}
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

BPDX=${BPDX:-64}
BPDY=${BPDY:-$((${BPDX}/2))} #${BPDY:-32}
BPDZ=${BPDZ:-$((${BPDX}/2))} #${BPDZ:-32}

NU=${NU:-0.00002840909091} # Re 550
BC=${BC:-freespace} # Re 550

FACTORY='IF3D_Sphere L=0.125 xpos=0.3 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1
'
# for accel and decel start and stop add accel=1 T=time_for_accel
# shift center to shed vortices immediately by ypos=0.250244140625 zpos=0.250244140625

OPTIONS=
OPTIONS+=" -bpdx ${BPDX} -bpdy ${BPDY} -bpdz ${BPDZ}"
OPTIONS+=" -dump2D 1 -dump3D 1 -tdump 0.5 -tend 8 -nslices 2 -slice1_direction 1 -slice2_direction 2 "
#OPTIONS+=" -BC_x dirichlet -BC_y dirichlet -BC_z dirichlet"
OPTIONS+=" -BC_x ${BC} -BC_y ${BC} -BC_z ${BC}"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1 -nu ${NU}"
