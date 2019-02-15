#!/bin/bash
NNODE=1
NNODEY=1
NNODE=$(($NNODEX * $NNODEY))

FACTORY='IF3D_Sphere L=0.125 xpos=0.3 xvel=0.125 bForcedInSimFrame=1 bFixFrameOfRef=1
'
# for accel and decel start and stop add accel=1 T=time_for_accel
# shift center to shed vortices immediately by ypos=0.250244140625 zpos=0.250244140625

OPTIONS=
#OPTIONS+=" -bpdx 64 -bpdy 32 -bpdz 32"
OPTIONS+=" -bpdx 32 -bpdy 16 -bpdz 16"
OPTIONS+=" -dump2D 0 -dump3D 1 -tdump 0.1 -tend 10"
OPTIONS+=" -BC_x dirichlet -BC_y dirichlet -BC_z dirichlet"
OPTIONS+=" -nprocsx ${NNODEX} -nprocsy ${NNODEY} -nprocsz 1"
OPTIONS+=" -CFL 0.1"
OPTIONS+=" -nu 0.00002840909091" # Re 550
